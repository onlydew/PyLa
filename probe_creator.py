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
DEFAULT_STRUT_MM = 0.0        # Default-Strebendicke (0 = keine Überlappung)

# Verzeichnisse
EXPORT_DIR_NAME = "export"
VARIANT_FOLDERS = [
    ("cell", "cell_variants", False),      # (Label, Ordner, allow_none?)
    ("pad",  "pad_variants", True),
    ("hook", "hook_variants", True),
    ("temp", "temp_variants", True),
    ("bed",  "bed_variants", True),
]


# Blender (für Boolean-Union Pads)
BLENDER_EXE = r"C:\Program Files\Blender Foundation\Blender 4.5\blender.exe"
BLENDER_BOOL_SCRIPT = "blender_boolean_operator.py"  # liegt im gleichen Ordner wie dieses Script

# Optionen
# interne Kontaktflächen auf Kachel-Ebenen entfernen (sinnvoll v.a. bei 0-Überlappung)
REMOVE_INTERNAL_INTERFACES = True
# echte doppelte Dreiecke entfernen (meist wenig Effekt, kann viel RAM ziehen)
DEDUPLICATE_TRIANGLES = True


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

FG_DIM = "\033[90m"

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

def find_cut_stl(path: str):
    """
    Zu einer gegebenen STL (…/foo.stl) die passende CUT-Version
    (…/foo_CUT.stl) suchen. Gibt Pfad zurück oder None.
    """
    root, ext = os.path.splitext(path)
    cut_path = root + "_CUT" + ext
    if os.path.exists(cut_path):
        return cut_path
    return None

def choose_variant_recursive(label: str, root_folder: str, allow_none: bool):
    """
    Interaktives Menü mit Pfeiltasten + Zahlen zur Auswahl einer STL-Datei.
    Navigiert rekursiv durch Unterordner.
    """
    import msvcrt

    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + f"  {label}-Variante auswählen" + RESET)
    print(SEPARATOR)

    if not os.path.isdir(root_folder):
        print(FG_RED + f"Ordner '{root_folder}' existiert nicht!" + RESET)
        if allow_none:
            print(FG_YELLOW + "Überspringe (keine Auswahl möglich)" + RESET)
            return None, None
        raise RuntimeError(f"Ordner '{root_folder}' fehlt – erforderlich.")

    root_abs = os.path.abspath(root_folder)
    current_dir = root_abs
    selected_index = 0

    while True:

        # -------- Inhalte auflisten --------
        entries = os.listdir(current_dir)
        subdirs = sorted([d for d in entries if os.path.isdir(os.path.join(current_dir, d))])
        stl_files = sorted([f for f in entries if f.lower().endswith(".stl")])

        items = [("dir", d) for d in subdirs] + [("file", f) for f in stl_files]

        # Extra-Optionen
        extra_items = []
        if current_dir != root_abs:
            extra_items.append(("back", ".."))
        if allow_none:
            extra_items.append(("none", "Keine Auswahl"))

        all_items = items + extra_items

        if not all_items:
            print(FG_YELLOW + "Keine Inhalte – gehe eine Ebene zurück." + RESET)
            current_dir = os.path.dirname(current_dir)
            selected_index = 0
            continue

        if selected_index >= len(all_items):
            selected_index = len(all_items) - 1

        # -------- Anzeige --------
        print("\033[2J\033[H", end="")  # Bildschirm löschen
        print(FG_MAGENTA + BOLD + f"{label}-Variante auswählen" + RESET)
        rel = os.path.relpath(current_dir, root_abs)
        print(FG_WHITE + "Ordner: " + FG_CYAN + rel + RESET)
        print()

        # normale Items
        for idx, (kind, name) in enumerate(items):
            prefix = FG_CYAN + ">" + RESET if idx == selected_index else " "
            if kind == "dir":
                print(f"{prefix} [Ordner] {name}")
            else:
                print(f"{prefix} [STL]    {name}")

        # extra items
        extra_start = len(items)
        for i, (kind, name) in enumerate(extra_items):
            idx = extra_start + i
            prefix = FG_CYAN + ">" + RESET if idx == selected_index else " "
            print(f"{prefix} {name}")

        print()
        print(FG_DIM + "↑/↓ bewegen | Enter auswählen | Zahl = direkte Auswahl | b = zurück" + RESET)

        # -------- Eingabe lesen --------
        ch = msvcrt.getch()

        # ENTER
        if ch in (b"\r", b"\n"):
            kind, name = all_items[selected_index]
            if kind == "dir":
                current_dir = os.path.join(current_dir, name)
                selected_index = 0
                continue
            elif kind == "file":
                full_path = os.path.join(current_dir, name)
                base_name = os.path.splitext(name)[0]
                return base_name, full_path
            elif kind == "back":
                current_dir = os.path.dirname(current_dir)
                selected_index = 0
                continue
            elif kind == "none":
                return None, None

        # Pfeiltasten
        elif ch in (b"\x00", b"\xe0"):
            ch2 = msvcrt.getch()
            if ch2 == b"H":  # ↑
                selected_index = (selected_index - 1) % len(all_items)
            elif ch2 == b"P":  # ↓
                selected_index = (selected_index + 1) % len(all_items)

        # Zahleneingabe
        elif ch.isdigit():
            digit = int(ch.decode())
            # 0 = keine Auswahl (wenn erlaubt)
            if digit == 0 and allow_none:
                return None, None
            # 1..n normale Items
            if 1 <= digit <= len(items):
                kind, name = items[digit - 1]
                if kind == "dir":
                    current_dir = os.path.join(current_dir, name)
                    selected_index = 0
                    continue
                else:
                    full_path = os.path.join(current_dir, name)
                    base_name = os.path.splitext(name)[0]
                    return base_name, full_path

        # „b“ für back
        elif ch.lower() == b"b":
            if current_dir == root_abs:
                continue
            current_dir = os.path.dirname(current_dir)
            selected_index = 0

        # andere Eingaben ignorieren
        else:
            pass

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

def scale_triangles(triangles, factor: float):
    """
    Skaliert alle Koordinaten eines Dreieck-Arrays um 'factor'.
    """
    if factor == 1.0:
        return triangles
    scaled = []
    for tri in triangles:
        scaled.append(tuple((x * factor, y * factor, z * factor) for (x, y, z) in tri))
    return scaled

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

    if sx == 0:
        sx = 1e-6
    if sy == 0:
        sy = 1e-6
    if sz == 0:
        sz = 1e-6

    normalized = [translate_triangle(tri, -minx, -miny, -minz)
                  for tri in triangles]

    return normalized, (sx, sy, sz)


def tile_cell(triangles, step, nx, ny, nz):
    """
    Kachelt die normalisierte Zelle nx * ny * nz mal in x,y,z-Richtung.

    step: (step_x, step_y, step_z) Abstand zwischen Zellursprüngen
          -> bei step_x < Zellkantenlänge entsteht eine Überlappung.
    """
    step_x, step_y, step_z = step
    tiled = []

    total_cells = nx * ny * nz
    iterator = tqdm(range(total_cells), desc="Tiling cells", unit="cell")

    for idx in iterator:
        ix = idx // (ny * nz)
        iy = (idx // nz) % ny
        iz = idx % nz

        dx = ix * step_x
        dy = iy * step_y
        dz = iz * step_z

        for tri in triangles:
            tiled.append(translate_triangle(tri, dx, dy, dz))

    return tiled


def remove_internal_interface_faces(triangles, step, nx, ny, nz, tol=1e-5):
    """
    Entfernt Dreiecke, die exakt auf internen Kachel-Ebenen liegen:

      x = k*step_x (1..nx-1),
      y = k*step_y (1..ny-1),
      z = k*step_z (1..nz-1)

    -> Das sind Kontaktflächen zwischen Zellen für den Fall
       Kante-an-Kante (kein Überlapp). Bei Überlappung lassen wir
       diese Funktion normalerweise aus.
    """
    print()
    print(FG_WHITE + "Entferne interne Kontaktflächen zwischen Zellen..." + RESET)

    step_x, step_y, step_z = step

    def on_internal_plane(coords, step_val, count):
        if count <= 1:
            return False
        k_float = coords[0] / step_val
        k = round(k_float)
        if k <= 0 or k >= count:
            return False
        plane = k * step_val
        return all(abs(c - plane) < tol for c in coords)

    kept = []
    removed = 0

    for tri in tqdm(triangles, desc="Removing internal faces", unit="tri"):
        xs = [v[0] for v in tri]
        ys = [v[1] for v in tri]
        zs = [v[2] for v in tri]

        remove = False

        if on_internal_plane(xs, step_x, nx):
            remove = True
        elif on_internal_plane(ys, step_y, ny):
            remove = True
        elif on_internal_plane(zs, step_z, nz):
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
        and "_cut.stl" not in f.lower()   # <-- CUT-Dateien ausblenden
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

def run_blender_boolean(op, mesh_a, mesh_b, output_stl):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, "blender_boolean_operator.py")

    if not os.path.exists(BLENDER_EXE):
        print(FG_RED + "Blender-Exe nicht gefunden." + RESET)
        return False
    if not os.path.exists(script_path):
        print(FG_RED + f"Blender-Script '{script_path}' nicht gefunden." + RESET)
        return False

    cmd = [
        BLENDER_EXE,
        "-b",
        "-P", script_path,
        "--",
        op,                 # 'UNION' oder 'DIFFERENCE'
        os.path.abspath(mesh_a),
        os.path.abspath(mesh_b),
        os.path.abspath(output_stl),
    ]

    print(FG_CYAN + f"Starte Blender Boolean ({op})..." + RESET)

    spinner = itertools.cycle("-\\|/")
    start_time = time.time()

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL,
    )

    # kleiner Spinner
    while True:
        ret = proc.poll()
        if ret is not None:
            break
        elapsed = time.time() - start_time
        print(
            f"\r{FG_CYAN}Boolean {op} läuft... {next(spinner)}  {elapsed:5.1f}s{RESET}",
            end="",
            flush=True
        )
        time.sleep(0.2)

    print()  # Zeilenumbruch

    out = proc.stdout.read().decode("utf-8", errors="ignore") if proc.stdout else ""
    log_path = os.path.join(os.path.dirname(output_stl), f"blender_{op.lower()}_log.txt")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(out)

    if proc.returncode != 0:
        print(FG_RED + f"Blender Boolean {op} fehlgeschlagen (Returncode {proc.returncode})." + RESET)
        print(FG_YELLOW + "Letzte Zeilen aus dem Blender-Log:" + RESET)
        print("\n".join(out.splitlines()[-20:]))
        return False

    print(FG_GREEN + f"Boolean {op} erfolgreich." + RESET)
    print(FG_WHITE + f"  Ausgabe: {FG_CYAN}{output_stl}{RESET}")
    print(FG_WHITE + f"  Log:     {FG_CYAN}{log_path}{RESET}")
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


def ask_float(prompt: str, default: float, min_value: float = 0.0) -> float:
    while True:
        s = input(f"{prompt} [{default}]: ").strip().replace(",", ".")
        if s == "":
            return default
        try:
            v = float(s)
            if v < min_value:
                print(FG_YELLOW + f"Wert muss >= {min_value} sein." + RESET)
                continue
            return v
        except ValueError:
            print(FG_YELLOW + "Bitte eine Zahl eingeben (z.B. 0.12)." + RESET)


# =====================================================
# Main
# =====================================================

def main():
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Lattice-Generator (direktes Tiling)" + RESET)
    print(SEPARATOR)

    # ------------------------------------------------------------------
    # 1) Varianten auswählen (Cell, Pad, Hook, Temp, Bed)
    # ------------------------------------------------------------------
    chosen_variants = {}       # label -> (name, path)
    variant_name_parts = []    # für den finalen Dateinamen

    for label, folder, allow_none in VARIANT_FOLDERS:
        nice_label = label.capitalize()
        name, path = choose_variant_recursive(nice_label, folder, allow_none)
        chosen_variants[label] = (name, path)

        if label == "cell":
            if name is None or path is None:
                raise RuntimeError("Es muss eine Cell-Variante gewählt werden.")
            variant_name_parts.append(name)
        else:
            if name is None:
                variant_name_parts.append(f"no{nice_label}")
            else:
                variant_name_parts.append(name)

    # Kurzreferenzen
    cell_name, cell_stl = chosen_variants["cell"]
    pad_name,  pad_stl  = chosen_variants["pad"]
    hook_name, hook_stl = chosen_variants["hook"]
    temp_name, temp_stl = chosen_variants["temp"]
    bed_name,  bed_stl  = chosen_variants["bed"]

    print()
    print(FG_WHITE + "Verwende Cell-Variante: " + FG_CYAN + cell_name + RESET)
    print(FG_WHITE + "Cell-STL:              " + FG_CYAN + cell_stl + RESET)

    # Hinweis zu den optionalen Varianten
    print()
    print(FG_WHITE + "Pad-Variante:   " + FG_CYAN + (pad_name  or "keine") + RESET)
    print(FG_WHITE + "Hook-Variante:  " + FG_CYAN + (hook_name or "keine") + RESET)
    print(FG_WHITE + "Temp-Variante:  " + FG_CYAN + (temp_name or "keine") + RESET)
    print(FG_WHITE + "Bed-Variante:   " + FG_CYAN + (bed_name  or "keine") + RESET)

    if not os.path.exists(cell_stl):
        raise FileNotFoundError(f"Eingabedatei '{cell_stl}' nicht gefunden.")

    # ------------------------------------------------------------------
    # 2) Anzahl Zellen pro Richtung
    # ------------------------------------------------------------------
    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Anzahl der Zellen pro Richtung" + RESET)
    print(SEPARATOR)

    nx = ask_int(FG_WHITE + "Anzahl Zellen in X" + RESET, DEFAULT_N_CELLS)
    ny = ask_int(FG_WHITE + "Anzahl Zellen in Y" + RESET, DEFAULT_N_CELLS)
    nz = ask_int(FG_WHITE + "Anzahl Zellen in Z" + RESET, DEFAULT_N_CELLS)

    # ------------------------------------------------------------------
    # 3) Export-Ordner + Basisname
    # ------------------------------------------------------------------
    export_dir = os.path.join(os.getcwd(), EXPORT_DIR_NAME)
    os.makedirs(export_dir, exist_ok=True)

    variant_suffix = "_".join(variant_name_parts)
    base_name = variant_suffix
    final_output_stl = os.path.join(export_dir, base_name + ".stl")

    # ------------------------------------------------------------------
    # 4) Einzelzell-STL lesen
    # ------------------------------------------------------------------
    print()
    print(FG_WHITE + f"Lese Cell-STL: {FG_CYAN}{cell_stl}{RESET}")
    base_triangles = read_stl_auto(cell_stl)
    base_tri_count = len(base_triangles)
    print(FG_WHITE + "  -> Dreiecke geladen: " + FG_CYAN + f"{base_tri_count:,}" + RESET)

    print_overview(base_tri_count, nx, ny, nz)

    print(FG_WHITE + "Normalisiere Zelle auf Ursprung und ermittle Zellgröße..." + RESET)
    norm_triangles, cell_size = normalize_to_origin(base_triangles)
    sx, sy, sz = cell_size
    print(FG_WHITE + f"  Zellgröße (Bounding Box): "
          + FG_CYAN + f"sx={sx:.3f} mm, sy={sy:.3f} mm, sz={sz:.3f} mm" + RESET)

    # ------------------------------------------------------------------
    # 5) Strebendicke + STL-Skalierung -> Schrittweiten
    # ------------------------------------------------------------------
    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + "  Strebendicke & STL-Skalierung" + RESET)
    print(SEPARATOR)

    # reale Strebendicke in mm
    strut_mm_user = ask_float(
        FG_WHITE + "Strebendicke in mm:" + RESET,
        DEFAULT_STRUT_MM,
        min_value=0.0
    )

    # STL-Skalierungsfaktor (z.B. 10.0, wenn Fusion 10x größer exportiert)
    stl_scale = ask_float(
        FG_WHITE + "Skalierungsfaktor der STL:" + RESET,
        1.0,
        min_value=0.0001
    )

    # Strebendicke im STL-Maßstab
    strut_mm_stl = strut_mm_user * stl_scale

    if strut_mm_user <= 0.0:
        # keine Überlappung -> Schritt = Zellkantenlänge
        step_x, step_y, step_z = sx, sy, sz
        print(FG_WHITE + "  Kante an Kante (keine Überlappung)." + RESET)
        print(FG_WHITE + "  Kachel-Schrittweiten: " +
              FG_CYAN + f"step_x={step_x:.3f}, step_y={step_y:.3f}, step_z={step_z:.3f} mm (STL-Raum)" + RESET)
    else:
        # Schritt = Zellkante - Strebendicke_im_STL
        step_x = max(sx - strut_mm_stl, sx * 0.1)
        step_y = max(sy - strut_mm_stl, sy * 0.1)
        step_z = max(sz - strut_mm_stl, sz * 0.1)

        if step_x != sx - strut_mm_stl or step_y != sy - strut_mm_stl or step_z != sz - strut_mm_stl:
            print(FG_YELLOW +
                  "Hinweis: Strebendicke*Skalierung war größer als 90% der Zellkante – "
                  "Schrittweite wurde zur Sicherheit begrenzt." + RESET)

        print(FG_WHITE + "  Eingestellte Strebendicke (real): " +
              FG_CYAN + f"{strut_mm_user:.4f} mm" + RESET)
        print(FG_WHITE + "  Strebendicke im STL-Maßstab:      " +
              FG_CYAN + f"{strut_mm_stl:.4f} mm" + RESET)
        print(FG_WHITE + "  Effektive Schrittweiten (Zellkante - Strebendicke_im_STL):" + RESET)
        print(FG_CYAN + f"    step_x={step_x:.4f}, step_y={step_y:.4f}, step_z={step_z:.4f} mm (STL-Raum)" + RESET)

    step = (step_x, step_y, step_z)

    # ------------------------------------------------------------------
    # 6) Tiling
    # ------------------------------------------------------------------
    print()
    print(FG_WHITE + f"Erzeuge Gitter: {FG_CYAN}{nx} x {ny} x {nz}" + RESET)
    lattice_triangles = tile_cell(norm_triangles, step, nx, ny, nz)
    print(FG_WHITE + "  -> Dreiecke nach Tiling: " +
          FG_CYAN + f"{len(lattice_triangles):,}" + RESET)

    # interne Kontaktflächen entfernen (nur sinnvoll bei strut_mm_user == 0)
    if REMOVE_INTERNAL_INTERFACES and strut_mm_user == 0.0:
        lattice_triangles = remove_internal_interface_faces(
            lattice_triangles, step, nx, ny, nz, tol=1e-5
        )

    # Duplikate entfernen (optional)
    if DEDUPLICATE_TRIANGLES:
        lattice_triangles = deduplicate_triangles(lattice_triangles, tol=1e-6)

     # ------------------------------------------------------------------
    # 7) Export / Boolean-Pipeline mit CUT-Unterstützung
    # ------------------------------------------------------------------
    success = False
    reported_path = None
    temp_files = []

    # Fall: es gibt überhaupt keine Zusatzgeometrie -> nur Lattice schreiben
    if all(stl is None for (_, stl) in [
        chosen_variants["pad"],
        chosen_variants["hook"],
        chosen_variants["temp"],
        chosen_variants["bed"],
    ]):
        print()
        print(FG_WHITE + f"Schreibe finales Binary-STL (nur Lattice): {FG_CYAN}{final_output_stl}{RESET}")
        # ggf. noch skalieren, falls du stl_scale verwendest:
        # scaled = scale_triangles(lattice_triangles, 1.0 / stl_scale)
        # write_binary_stl(scaled, final_output_stl, base_name, show_progress=True)
        write_binary_stl(lattice_triangles, final_output_stl, base_name, show_progress=True)
        print(FG_GREEN + "STL-Export fertig." + RESET)
        success = True
        reported_path = final_output_stl

    else:
        # 1) Lattice-only STL schreiben (Zwischenschritt)
        lattice_only_path = os.path.join(export_dir, base_name + "_step0_lattice.stl")
        print()
        print(FG_WHITE + f"Schreibe Lattice-only STL für Boolean-Pipeline: {FG_CYAN}{lattice_only_path}{RESET}")
        write_binary_stl(lattice_triangles, lattice_only_path, base_name, show_progress=True)
        temp_files.append(lattice_only_path)

        current = lattice_only_path

        # Reihenfolge der Varianten, für die Booleans laufen sollen
        pipeline = [
            ("pad",  pad_stl),
            ("hook", hook_stl),
            ("temp", temp_stl),
            ("bed",  bed_stl),
        ]

        for label, stl_path in pipeline:
            if stl_path is None:
                continue

            nice = label.capitalize()
            print()
            print(FG_CYAN + f"=== Boolean mit {nice} ===" + RESET)

            # 1) CUT-Datei suchen (…_CUT.stl)
            cut_path = find_cut_stl(stl_path)
            if cut_path:
                diff_out = os.path.join(export_dir, f"{base_name}_diff_{label}.stl")
                print(FG_WHITE + f"  CUT gefunden: {FG_CYAN}{cut_path}{RESET}")
                print(FG_WHITE + f"  Führe DIFFERENCE aus: current - CUT -> {diff_out}" + RESET)
                ok = run_blender_boolean("DIFFERENCE", current, cut_path, diff_out)
                if not ok:
                    raise RuntimeError(f"Boolean DIFFERENCE mit {label} fehlgeschlagen.")
                temp_files.append(diff_out)
                current = diff_out

            # 2) UNION mit der eigentlichen STL
            union_out = os.path.join(export_dir, f"{base_name}_union_{label}.stl")
            print(FG_WHITE + f"  Führe UNION aus: current ∪ {stl_path} -> {union_out}" + RESET)
            ok = run_blender_boolean("UNION", current, stl_path, union_out)
            if not ok:
                raise RuntimeError(f"Boolean UNION mit {label} fehlgeschlagen.")
            temp_files.append(union_out)
            current = union_out

        # Am Ende liegt das letzte Boolean-Ergebnis in `current`
        final_boolean_stl = current

        print()
        print(FG_WHITE + "Übernehme finales Boolean-Ergebnis als Output-STL..." + RESET)
        # Wenn du noch skalieren willst, hier STL einlesen + skaliert schreiben,
        # ansonsten einfach Datei kopieren:
        if final_boolean_stl != final_output_stl:
            # einfache Kopie
            import shutil
            shutil.copy2(final_boolean_stl, final_output_stl)

        success = True
        reported_path = final_output_stl

    # ------------------------------------------------------------------
    # 8) temporäre Dateien aufräumen
    # ------------------------------------------------------------------
    print()
    print(FG_WHITE + "Bereinige temporäre Dateien..." + RESET)
    for path in set(temp_files):
        try:
            if os.path.exists(path):
                os.remove(path)
                print(FG_DIM + f"  gelöscht: {path}" + RESET)
        except OSError as e:
            print(FG_YELLOW + f"  konnte {path} nicht löschen: {e}" + RESET)

    # ------------------------------------------------------------------
    # 9) Abschluss
    # ------------------------------------------------------------------
    print()
    print(SEPARATOR)
    print(FG_GREEN + BOLD + "  Alles fertig." + RESET)
    print(SEPARATOR)
    if os.path.exists(final_output_stl):
        print(FG_WHITE + "  Ausgabe-STL: " + FG_CYAN + final_output_stl + RESET)
    else:
        print(FG_WHITE + "  Finale STL konnte NICHT erzeugt werden." + RESET)
    print(SEPARATOR)

if __name__ == "__main__":
    main()
