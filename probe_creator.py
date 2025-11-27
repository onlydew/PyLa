import os
import struct
import subprocess
import time
import itertools
from tqdm import tqdm

# =====================================================
# Konfiguration
# =====================================================

# Defaultwerte für interaktive Eingaben
DEFAULT_N_CELLS = 20          # Default für Nx, Ny, Nz

# Verzeichnisse
CELL_TYPES_DIR = "cell_types"
PAD_TYPES_DIR = "pad_types"
EXPORT_DIR_NAME = "export"

# Blender (für Boolean-Union Pads)
BLENDER_EXE = r"C:\Program Files\Blender Foundation\Blender 4.4\blender.exe"
BLENDER_BOOL_SCRIPT = "blender_boolean_union.py"  # liegt im gleichen Ordner wie dieses Script

# Optionen
REMOVE_INTERNAL_INTERFACES = True   # interne Kontaktflächen zwischen Zellen entfernen
DEDUPLICATE_TRIANGLES = False       # echte doppelte Dreiecke entfernen (meist wenig Effekt, kann viel RAM ziehen)


# =====================================================
# ANSI-Farben für hübschere Ausgabe
# =====================================================

RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"

FG_RED = "\033[31m"
FG_GREEN = "\033[32m"
FG_YELLOW = "\033[33m"
FG_BLUE = "\033[34m"
FG_MAGENTA = "\033[35m"
FG_CYAN = "\033[36m"
FG_WHITE = "\033[97m"

SEPARATOR = FG_BLUE + "-" * 72 + RESET


# =====================================================
# STL-Hilfsfunktionen (ASCII + Binary)
# =====================================================

def read_ascii_stl(path: str):
    triangles = []
    current_triangle = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if parts[0].lower() == "vertex" and len(parts) >= 4:
                x, y, z = map(float, parts[1:4])
                current_triangle.append((x, y, z))

                if len(current_triangle) == 3:
                    triangles.append(tuple(current_triangle))
                    current_triangle = []

    if not triangles:
        raise ValueError(
            f"Keine Dreiecke in {path} gefunden – ist es wirklich ASCII STL?"
        )
    return triangles


def read_binary_stl(path: str):
    triangles = []
    with open(path, "rb") as f:
        header = f.read(80)
        n_bytes = f.read(4)
        if len(n_bytes) != 4:
            raise ValueError(f"Datei {path} zu kurz für gültige Binary-STL.")
        (num_triangles,) = struct.unpack("<I", n_bytes)

        for _ in range(num_triangles):
            data = f.read(50)
            if len(data) != 50:
                raise ValueError(
                    f"Datei {path}: Unerwartetes Ende beim Lesen der Dreiecke."
                )
            unpacked = struct.unpack("<12fH", data)
            v1 = (unpacked[3], unpacked[4], unpacked[5])
            v2 = (unpacked[6], unpacked[7], unpacked[8])
            v3 = (unpacked[9], unpacked[10], unpacked[11])
            triangles.append((v1, v2, v3))

    if not triangles:
        raise ValueError(f"Keine Dreiecke in {path} gefunden – Binary-STL leer?")
    return triangles


def read_stl_auto(path: str):
    with open(path, "rb") as f:
        start = f.read(512)
    text = start.decode("ascii", errors="ignore").lower()

    if text.lstrip().startswith("solid") and "facet" in text:
        print(FG_CYAN + "  Erkanntes Format: ASCII STL" + RESET)
        return read_ascii_stl(path)
    else:
        print(FG_CYAN + "  Erkanntes Format: Binary STL" + RESET)
        return read_binary_stl(path)


def write_binary_stl(triangles, path, solid_name="mesh", show_progress=True):
    with open(path, "wb") as f:
        header = (solid_name[:79]).ljust(80, " ").encode("ascii")
        f.write(header)

        f.write(struct.pack("<I", len(triangles)))

        iterable = triangles
        if show_progress:
            iterable = tqdm(triangles, desc="Writing Binary STL", unit="tri")

        for tri in iterable:
            data = struct.pack(
                "<12fH",
                0.0, 0.0, 0.0,                  # Normal (Dummy)
                tri[0][0], tri[0][1], tri[0][2],
                tri[1][0], tri[1][1], tri[1][2],
                tri[2][0], tri[2][1], tri[2][2],
                0                               # Attribute Byte Count
            )
            f.write(data)


# =====================================================
# Geometrie-Hilfsfunktionen
# =====================================================

def compute_bounding_box(triangles):
    xs, ys, zs = [], [], []
    for tri in triangles:
        for (x, y, z) in tri:
            xs.append(x)
            ys.append(y)
            zs.append(z)
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def translate_triangle(tri, dx, dy, dz):
    return tuple((x + dx, y + dy, z + dz) for (x, y, z) in tri)


def normalize_to_origin(triangles):
    """
    Verschiebt die Zelle so, dass ihre Bounding Box bei (0,0,0) beginnt.
    Gibt neue Dreiecke + Zellgröße (sx, sy, sz) zurück.
    """
    minx, maxx, miny, maxy, minz, maxz = compute_bounding_box(triangles)
    sx = maxx - minx
    sy = maxy - miny
    sz = maxz - minz

    if sx == 0: sx = 1e-6
    if sy == 0: sy = 1e-6
    if sz == 0: sz = 1e-6

    normalized = [translate_triangle(tri, -minx, -miny, -minz)
                  for tri in triangles]

    return normalized, (sx, sy, sz)


def tile_cell(triangles, cell_size, nx, ny, nz):
    """
    Kachelt die normalisierte Zelle nx * ny * nz mal in x,y,z-Richtung.
    """
    sx, sy, sz = cell_size
    tiled = []

    total_cells = nx * ny * nz
    iterator = tqdm(range(total_cells), desc="Tiling cells", unit="cell")

    for idx in iterator:
        ix = idx // (ny * nz)
        iy = (idx // nz) % ny
        iz = idx % nz

        dx = ix * sx
        dy = iy * sy
        dz = iz * sz

        for tri in triangles:
            tiled.append(translate_triangle(tri, dx, dy, dz))

    return tiled


def remove_internal_interface_faces(triangles, cell_size, nx, ny, nz, tol=1e-5):
    """
    Entfernt Dreiecke, die exakt auf internen Kachel-Ebenen liegen:
      x = k*sx (1..nx-1), y = k*sy (1..ny-1), z = k*sz (1..nz-1)
    -> Das sind die Kontaktflächen zwischen benachbarten Zellen.
    Außenflächen (x=0, x=nx*sx usw.) bleiben erhalten.
    """
    print()
    print(FG_WHITE + "Entferne interne Kontaktflächen zwischen Zellen..." + RESET)

    sx, sy, sz = cell_size

    def on_internal_plane(coords, step, count):
        k_float = coords[0] / step
        k = round(k_float)
        if k <= 0 or k >= count:
            return False
        plane = k * step
        return all(abs(c - plane) < tol for c in coords)

    kept = []
    removed = 0

    for tri in tqdm(triangles, desc="Removing internal faces", unit="tri"):
        xs = [v[0] for v in tri]
        ys = [v[1] for v in tri]
        zs = [v[2] for v in tri]

        remove = False

        if nx > 1 and on_internal_plane(xs, sx, nx):
            remove = True
        elif ny > 1 and on_internal_plane(ys, sy, ny):
            remove = True
        elif nz > 1 and on_internal_plane(zs, sz, nz):
            remove = True

        if remove:
            removed += 1
        else:
            kept.append(tri)

    print(FG_WHITE + "  Vorher:  " + FG_CYAN + f"{len(triangles):,}" + RESET)
    print(FG_WHITE + "  Nachher: " + FG_CYAN + f"{len(kept):,}" + RESET)
    print(FG_WHITE + "  Entfernt:" + FG_CYAN + f" {removed:,}" + RESET)
    return kept


def deduplicate_triangles(triangles, tol=1e-6):
    """
    Entfernt doppelte Dreiecke (gleiche 3 Punkte, Reihenfolge egal).
    Achtung: Kann sehr speicherintensiv werden bei vielen Dreiecken.
    """
    print()
    print(FG_WHITE + "Entferne doppelte Dreiecke..." + RESET)
    inv_tol = 1.0 / tol

    def quantize(p):
        return (
            int(round(p[0] * inv_tol)),
            int(round(p[1] * inv_tol)),
            int(round(p[2] * inv_tol)),
        )

    seen = set()
    unique = []

    for tri in tqdm(triangles, desc="Deduplicating", unit="tri"):
        qv = tuple(sorted(quantize(v) for v in tri))
        if qv not in seen:
            seen.add(qv)
            unique.append(tri)

    print(FG_WHITE + "  Vorher:  " + FG_CYAN + f"{len(triangles):,}" + RESET)
    print(FG_WHITE + "  Nachher: " + FG_CYAN + f"{len(unique):,}" + RESET)
    print(FG_WHITE + "  Entfernt:" + FG_CYAN + f" {len(triangles) - len(unique):,}" + RESET)
    return unique


# =====================================================
# Übersicht
# =====================================================

def print_overview(base_tri_count, nx, ny, nz):
    cells = nx * ny * nz
    est_tris = base_tri_count * cells
    est_bytes_binary = est_tris * 50  # 50 Bytes pro Dreieck in Binary STL
    est_mb = est_bytes_binary / (1024 ** 2)

    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Übersicht (grobe Abschätzung)" + RESET)
    print(SEPARATOR)
    print(f"{FG_WHITE}  Dreiecke pro Zelle: {FG_CYAN}{base_tri_count:,}{RESET}")
    print(f"{FG_WHITE}  Zellen im Würfel:   {FG_CYAN}{cells:,}{RESET}"
          f"{FG_WHITE}  ({nx} x {ny} x {nz}){RESET}")
    print(f"{FG_WHITE}  Max. Dreiecke roh:  {FG_CYAN}{est_tris:,}{RESET}")
    print(f"{FG_WHITE}  Rohes Binary-STL:   {FG_CYAN}~{est_mb:.2f} MB{RESET}")
    print(SEPARATOR)
    print()


# =====================================================
# Typen aus Ordnern lesen (Cells & Pads)
# =====================================================

def discover_stl_types_from_folder(folder: str):
    if not os.path.isdir(folder):
        return {}

    stl_files = sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith(".stl")
    )

    types = {}
    idx = 1
    for fname in stl_files:
        name = os.path.splitext(fname)[0]
        full_path = os.path.join(folder, fname)
        types[str(idx)] = (name, full_path)
        idx += 1

    return types


def choose_cell_type():
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Zelltyp auswählen (aus Ordner 'cell_types/')" + RESET)
    print(SEPARATOR)

    cell_types = discover_stl_types_from_folder(CELL_TYPES_DIR)
    if not cell_types:
        raise RuntimeError(
            f"Keine STL-Dateien in '{CELL_TYPES_DIR}' gefunden. "
            f"Jede .stl-Datei entspricht einem Zelltyp."
        )

    for key, (name, stl_path) in cell_types.items():
        print(f"  {FG_CYAN}{key}{RESET}: {FG_WHITE}{name}{RESET}   "
              f"{DIM}(Datei: {os.path.basename(stl_path)}){RESET}")

    print()

    while True:
        choice = input(FG_WHITE + "Deine Auswahl: " + RESET).strip()
        if choice in cell_types:
            return cell_types[choice]
        print(FG_YELLOW + "Ungültige Auswahl. Bitte erneut versuchen." + RESET)


def choose_pad_type():
    """
    Wählt ein Pad aus 'pad_types/' aus. Gibt (pad_name, pad_path) zurück oder (None, None),
    wenn kein Pad verwendet werden soll.
    """
    pad_types = discover_stl_types_from_folder(PAD_TYPES_DIR)

    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Pad-Typ auswählen (aus 'pad_types/')" + RESET)
    print(SEPARATOR)

    print(f"  {FG_CYAN}0{RESET}: {FG_WHITE}Kein Pad verwenden{RESET}")

    if not pad_types:
        print(FG_YELLOW + "  (Keine Pad-STLs gefunden – '0' wird automatisch gewählt.)" + RESET)
        return None, None

    for key, (name, stl_path) in pad_types.items():
        print(f"  {FG_CYAN}{key}{RESET}: {FG_WHITE}{name}{RESET}   "
              f"{DIM}(Datei: {os.path.basename(stl_path)}){RESET}")

    print()

    while True:
        choice = input(FG_WHITE + "Deine Auswahl: " + RESET).strip()
        if choice == "0":
            return None, None
        if choice in pad_types:
            return pad_types[choice]
        print(FG_YELLOW + "Ungültige Auswahl. Bitte erneut versuchen." + RESET)


# =====================================================
# Blender Headless: Boolean Union (Lattice ∪ Pads)
# =====================================================

def run_blender_boolean_union(lattice_stl: str, pads_stl: str, output_stl: str) -> bool:
    """
    Startet Blender im Hintergrund und führt die Boolean-Union-Operation aus:
      result = lattice ∪ pads

    Erwartet, dass 'blender_boolean_union.py' im selben Ordner
    wie dieses Script liegt.
    Schreibt die Blender-Konsole zusätzlich in eine Logdatei, um
    Fehler besser debuggen zu können.
    """
    # --- Pfade robust bestimmen ---
    script_dir = os.path.dirname(os.path.abspath(__file__))  # Ordner von probe_creator.py
    script_path = os.path.join(script_dir, BLENDER_BOOL_SCRIPT)

    if not os.path.exists(BLENDER_EXE):
        print(FG_YELLOW + "Blender-Exe nicht gefunden – Boolean-Union übersprungen." + RESET)
        print(FG_YELLOW + "Pfad in BLENDER_EXE anpassen." + RESET)
        return False

    if not os.path.exists(script_path):
        print(FG_RED + f"Blender-Script '{BLENDER_BOOL_SCRIPT}' nicht gefunden." + RESET)
        print(FG_RED + "Bitte im selben Ordner wie 'probe_creator.py' ablegen." + RESET)
        print(FG_RED + f"Gesuchter Pfad war: {script_path}" + RESET)
        return False

    cmd = [
        BLENDER_EXE,
        "-b",
        "-P", script_path,
        "--",
        os.path.abspath(lattice_stl),
        os.path.abspath(pads_stl),
        os.path.abspath(output_stl),
    ]

    print()
    print(FG_CYAN + "Starte Blender für Boolean-Union (Lattice ∪ Pads)..." + RESET)

    # Logfile im Export-Ordner
    export_dir = os.path.join(script_dir, EXPORT_DIR_NAME)
    os.makedirs(export_dir, exist_ok=True)
    log_path = os.path.join(export_dir, "blender_boolean_union.log")

    spinner = itertools.cycle("-\\|/")
    start_time = time.time()

    try:
        # Blender-Ausgabe in Logfile schreiben, damit wir sie bei Fehlern lesen können
        log_file = open(log_path, "wb")
        proc = subprocess.Popen(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        print(FG_RED + "Blender konnte nicht gestartet werden. Pfad prüfen." + RESET)
        return False
    except Exception as e:
        print(FG_RED + f"Fehler beim Starten von Blender: {e}" + RESET)
        return False

    # Spinner anzeigen, solange Blender läuft
    while True:
        ret = proc.poll()
        if ret is not None:
            break
        elapsed = time.time() - start_time
        print(
            f"\r{FG_CYAN}Boolean-Union läuft... {next(spinner)}  {elapsed:5.1f}s{RESET}",
            end="",
            flush=True
        )
        time.sleep(0.2)

    # Logfile schließen, bevor wir es evtl. lesen
    log_file.close()
    print()  # Zeilenumbruch nach dem Spinner

    if proc.returncode != 0:
        print(FG_RED + f"Blender Boolean-Union fehlgeschlagen (Returncode {proc.returncode})." + RESET)
        print(FG_YELLOW + "Blender-Log findest du hier:" + RESET)
        print(FG_CYAN + log_path + RESET)

        # Optional: letzte paar Zeilen direkt anzeigen
        try:
            with open(log_path, "r", encoding="utf-8", errors="ignore") as lf:
                lines = lf.readlines()
            tail = "".join(lines[-20:])  # letzte 20 Zeilen
            print(FG_YELLOW + "Letzte Zeilen aus dem Blender-Log:" + RESET)
            print(FG_WHITE + tail + RESET)
        except Exception as e:
            print(FG_YELLOW + f"(Konnte Log nicht lesen: {e})" + RESET)

        return False

    # Erfolgspfad
    if not os.path.exists(output_stl):
        print(FG_RED + "Blender lief durch, aber die Output-STL existiert nicht!" + RESET)
        print(FG_YELLOW + "Pfad, den wir erwartet haben:" + RESET)
        print(FG_CYAN + output_stl + RESET)
        print(FG_YELLOW + "Blender-Log zur Analyse:" + RESET)
        print(FG_CYAN + log_path + RESET)
        return False

    elapsed = time.time() - start_time
    size_bytes = os.path.getsize(output_stl)
    print(FG_GREEN + f"Boolean-Union fertig in {elapsed:5.1f}s." + RESET)
    print(FG_GREEN + f"Boolean-Ergebnis geschrieben nach: {output_stl}" + RESET)
    print(FG_WHITE + f"  Dateigröße: {FG_CYAN}{size_bytes/1_000_000:.2f} MB" + RESET)
    print(FG_WHITE + f"  Blender-Log: {FG_CYAN}{log_path}{RESET}")

    return True

# =====================================================
# Interaktive Eingaben
# =====================================================

def ask_int(prompt: str, default: int, min_value: int = 1) -> int:
    while True:
        s = input(f"{prompt} [{default}]: ").strip()
        if s == "":
            return default
        try:
            v = int(s)
            if v < min_value:
                print(FG_YELLOW + f"Wert muss >= {min_value} sein." + RESET)
                continue
            return v
        except ValueError:
            print(FG_YELLOW + "Bitte eine ganzzahlige Zahl eingeben." + RESET)


# =====================================================
# Main
# =====================================================

def main():
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Lattice-Generator" + RESET)
    print(SEPARATOR)

    # --- Zelltyp auswählen ---
    cell_name, input_stl = choose_cell_type()
    print()
    print(FG_WHITE + "Verwende Zelltyp: " + FG_CYAN + cell_name + RESET)
    print(FG_WHITE + "Einzelzell-STL:   " + FG_CYAN + input_stl + RESET)

    if not os.path.exists(input_stl):
        raise FileNotFoundError(f"Eingabedatei '{input_stl}' nicht gefunden.")

    # --- Anzahl Zellen pro Richtung abfragen ---
    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Anzahl der Zellen pro Richtung" + RESET)
    print(SEPARATOR)

    nx = ask_int(FG_WHITE + "Anzahl Zellen in X" + RESET, DEFAULT_N_CELLS)
    ny = ask_int(FG_WHITE + "Anzahl Zellen in Y" + RESET, DEFAULT_N_CELLS)
    nz = ask_int(FG_WHITE + "Anzahl Zellen in Z" + RESET, DEFAULT_N_CELLS)

    # --- Pad auswählen (optional) ---
    pad_name, pad_stl = choose_pad_type()
    if pad_name is not None:
        print()
        print(FG_WHITE + "Verwende Pad-Typ: " + FG_CYAN + pad_name + RESET)
        print(FG_WHITE + "Pad-STL:         " + FG_CYAN + pad_stl + RESET)

    # Export-Ordner erstellen
    export_dir = os.path.join(os.getcwd(), EXPORT_DIR_NAME)
    os.makedirs(export_dir, exist_ok=True)

    # Basisname für Dateien
    base_name = f"{cell_name}_{nx}x{ny}x{nz}"
    if pad_name is not None:
        safe_pad_name = pad_name.replace(" ", "-")
        base_name += f"_{safe_pad_name}"

    # Finale Ausgabe
    final_output_stl = os.path.join(export_dir, base_name + ".stl")

    # --- Einzelzell-STL laden ---
    print()
    print(FG_WHITE + f"Lese Einzelzell-STL: {FG_CYAN}{input_stl}{RESET}")
    base_triangles = read_stl_auto(input_stl)
    base_tri_count = len(base_triangles)
    print(FG_WHITE + "  -> Dreiecke geladen: " + FG_CYAN + f"{base_tri_count:,}" + RESET)

    print_overview(base_tri_count, nx, ny, nz)

    print(FG_WHITE + "Normalisiere Zelle auf Ursprung und ermittle Zellgröße..." + RESET)
    norm_triangles, cell_size = normalize_to_origin(base_triangles)
    sx, sy, sz = cell_size
    print(FG_WHITE + f"  Zellgröße: {FG_CYAN}sx={sx:.3f} mm, sy={sy:.3f} mm, sz={sz:.3f} mm" + RESET)

    # --- Tiling ---
    print()
    print(FG_WHITE + f"Erzeuge Gitter: {FG_CYAN}{nx} x {ny} x {nz}" + RESET)
    lattice_triangles = tile_cell(norm_triangles, cell_size, nx, ny, nz)
    print(FG_WHITE + "  -> Dreiecke nach Tiling: " +
          FG_CYAN + f"{len(lattice_triangles):,}" + RESET)

    # --- interne Kontaktflächen entfernen (nur im Lattice) ---
    if REMOVE_INTERNAL_INTERFACES:
        lattice_triangles = remove_internal_interface_faces(
            lattice_triangles, cell_size, nx, ny, nz, tol=1e-5
        )

    # --- Duplikate entfernen (optional) ---
    if DEDUPLICATE_TRIANGLES:
        lattice_triangles = deduplicate_triangles(lattice_triangles, tol=1e-6)

    # Erfolg / Pfad für Abschlussausgabe
    success = False
    reported_path = None

    # --- Fall 1: keine Pads -> direkt final schreiben ---
    if pad_name is None or pad_stl is None:
        print()
        print(FG_WHITE + f"Schreibe finales Binary-STL (ohne Pads): {FG_CYAN}{final_output_stl}{RESET}")
        write_binary_stl(lattice_triangles, final_output_stl, base_name, show_progress=True)
        print(FG_GREEN + "STL-Export fertig." + RESET)

        success = True
        reported_path = final_output_stl

    # --- Fall 2: mit Pads -> erst Lattice-only schreiben, dann Boolean in Blender ---
    else:
        lattice_only_path = os.path.join(export_dir, base_name + "_lattice_only_tmp.stl")

        print()
        print(FG_WHITE + f"Schreibe Lattice-only STL für Boolean-Union: {FG_CYAN}{lattice_only_path}{RESET}")
        write_binary_stl(lattice_triangles, lattice_only_path, base_name, show_progress=True)

        # Boolean-Union: (Lattice-only) ∪ (Pad-STL) -> final_output_stl
        ok = run_blender_boolean_union(lattice_only_path, pad_stl, final_output_stl)

        if ok:
            # Boolean erfolgreich -> Temp kann weg
            try:
                os.remove(lattice_only_path)
                print(FG_WHITE + "Temporäre Lattice-only STL gelöscht." + RESET)
            except OSError as e:
                print(FG_YELLOW + f"Konnte temporäre Lattice-only STL nicht löschen: {e}" + RESET)

            print(FG_GREEN + "STL-Export (Lattice ∪ Pads) fertig." + RESET)
            success = True
            reported_path = final_output_stl
        else:
            # Boolean fehlgeschlagen -> Lattice-only als Fallback behalten
            print(FG_YELLOW + "Boolean-Union fehlgeschlagen – Lattice-only STL bleibt als Fallback erhalten." + RESET)
            print(FG_YELLOW + "Du findest das Lattice ohne Pads unter:" + RESET)
            print(FG_CYAN + lattice_only_path + RESET)
            success = False
            reported_path = lattice_only_path

    # --- Abschluss ---
    print()
    print(SEPARATOR)
    print(FG_GREEN + BOLD + "  Alles fertig." + RESET)
    print(SEPARATOR)
    if reported_path is not None and os.path.exists(reported_path):
        print(FG_WHITE + "  Ausgabe-STL: " + FG_CYAN + reported_path + RESET)
    else:
        print(FG_WHITE + "  Finale STL konnte NICHT erzeugt werden." + RESET)
    print(SEPARATOR)

if __name__ == "__main__":
    main()
