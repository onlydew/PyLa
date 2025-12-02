import os
import time
import itertools

from tqdm import tqdm

import cadquery as cq
from cadquery import exporters, importers

# =====================================================
# Konfiguration
# =====================================================

# Defaultwerte für interaktive Eingaben
DEFAULT_N_CELLS = 20      # Default für Nx, Ny, Nz
DEFAULT_STRUT_MM = 0.0    # Default-Strebendicke (0 = keine Überlappung)

# Verzeichnisse
EXPORT_DIR_NAME = "export"

VARIANT_FOLDERS = [
    ("cell", "cell_variants", False),  # (Label, Ordner, allow_none?)
    ("pad",  "pad_variants",  True),
    ("hook", "hook_variants", True),
    ("temp", "temp_variants", True),
    ("bed",  "bed_variants",  True),
]

# Optionen (einige STL-spezifische Optionen entfallen)
# Hier könntest du später noch CAD-spezifische Optimierungen ergänzen

# =====================================================
# ANSI-Farben für hübschere Ausgabe
# =====================================================

RESET    = "\033[0m"
BOLD     = "\033[1m"
DIM      = "\033[2m"

FG_RED    = "\033[31m"
FG_GREEN  = "\033[32m"
FG_YELLOW = "\033[33m"
FG_BLUE   = "\033[34m"
FG_MAGENTA= "\033[35m"
FG_CYAN   = "\033[36m"
FG_WHITE  = "\033[97m"
FG_DIM    = "\033[90m"

SEPARATOR = FG_BLUE + "-" * 72 + RESET

# =====================================================
# Datei-Auswahl (wie früher, aber für STEP statt STL)
# =====================================================

def choose_variant_recursive(label: str, root_folder: str, allow_none: bool):
    """
    Interaktives Menü mit Pfeiltasten + Zahlen zur Auswahl einer STEP-Datei.
    Navigiert rekursiv durch Unterordner.

    Erwartet .step / .stp Dateien.
    """
    import msvcrt

    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + f" {label}-Variante auswählen" + RESET)
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
        subdirs = sorted(
            [d for d in entries if os.path.isdir(os.path.join(current_dir, d))]
        )
        step_files = sorted(
            [
                f for f in entries
                if f.lower().endswith(".step") or f.lower().endswith(".stp")
            ]
        )

        items = [("dir", d) for d in subdirs] + [("file", f) for f in step_files]

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
                print(f"{prefix} [STEP]   {name}")

        # extra items
        extra_start = len(items)
        for i, (kind, name) in enumerate(extra_items):
            idx = extra_start + i
            prefix = FG_CYAN + ">" + RESET if idx == selected_index else " "
            print(f"{prefix} {name}")

        print()
        print(
            FG_DIM
            + "↑/↓ bewegen | Enter auswählen | Zahl = direkte Auswahl | b = zurück"
            + RESET
        )

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
            if ch2 == b"H":   # ↑
                selected_index = (selected_index - 1) % len(all_items)
            elif ch2 == b"P": # ↓
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
# CAD-Hilfsfunktionen (CadQuery)
# =====================================================

def load_step_shape(path: str):
    """
    Lädt eine STEP-Datei mit CadQuery und gibt ein Shape zurück.
    """
    shape_or_wp = importers.importStep(path)
    # CadQuery kann entweder ein Shape oder eine Workplane liefern
    if isinstance(shape_or_wp, cq.Workplane):
        shape = shape_or_wp.val()
    else:
        shape = shape_or_wp
    return shape

def normalize_shape_to_origin(shape: cq.Shape):
    """
    Verschiebt die Zelle so, dass ihre Bounding Box bei (0,0,0) beginnt.
    Gibt neues Shape + Zellgröße (sx, sy, sz) zurück.
    """
    bb = shape.BoundingBox()
    tx = -bb.xmin
    ty = -bb.ymin
    tz = -bb.zmin

    loc = cq.Location(cq.Vector(tx, ty, tz))
    norm_shape = shape.located(loc)

    sx = bb.xlen
    sy = bb.ylen
    sz = bb.zlen

    # kleine Sicherheit gegen 0-Kantenlänge
    eps = 1e-6
    sx = sx if sx != 0 else eps
    sy = sy if sy != 0 else eps
    sz = sz if sz != 0 else eps

    return norm_shape, (sx, sy, sz)

def build_lattice_cell_tiling(cell_shape: cq.Shape, step, nx, ny, nz):
    """
    Baut ein Nx * Ny * Nz Lattice aus einer normalisierten Zelle.
    step: (step_x, step_y, step_z) Abstand zwischen Zellursprüngen.
    Gibt ein CadQuery-Compound (Shape) zurück.
    """
    step_x, step_y, step_z = step

    asm = cq.Assembly()
    total_cells = nx * ny * nz
    iterator = tqdm(range(total_cells), desc="Tiling cells", unit="cell")

    for idx in iterator:
        ix = idx // (ny * nz)
        iy = (idx // nz) % ny
        iz = idx % nz

        dx = ix * step_x
        dy = iy * step_y
        dz = iz * step_z

        loc = cq.Location(cq.Vector(dx, dy, dz))
        # Wir fügen immer dieselbe Zelle an unterschiedlichen Positionen ein
        asm.add(cell_shape, loc=loc)

    # Assembly zu einem Compound zusammenfassen
    compound = asm.toCompound()
    return compound

def union_optional_variant(main_shape: cq.Shape, variant_step: str, label: str):
    """
    Lädt eine optionale STEP-Variante (pad/hook/temp/bed),
    und verschmilzt sie per Boolean-Union mit dem Haupt-Shape.
    """
    if variant_step is None:
        return main_shape

    print()
    print(FG_CYAN + f"=== Boolean-Union mit {label} ===" + RESET)
    print(FG_WHITE + f" Lade STEP: {FG_CYAN}{variant_step}{RESET}")

    var_shape = load_step_shape(variant_step)

    # Hinweis: Es wird angenommen, dass die Variante bereits in den
    # gleichen Koordinaten wie das Lattice modelliert wurde.
    start_time = time.time()
    fused = main_shape.fuse(var_shape)
    elapsed = time.time() - start_time

    print(FG_GREEN + f" Union mit {label} fertig in {elapsed:.1f}s." + RESET)
    return fused

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
# Übersicht (Abschätzung – hier nur Zellanzahl)
# =====================================================

def print_overview(nx, ny, nz):
    cells = nx * ny * nz
    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + " Übersicht (grobe Abschätzung)" + RESET)
    print(SEPARATOR)
    print(f"{FG_WHITE} Zellen im Würfel: {FG_CYAN}{cells:,}{RESET}"
          f"{FG_WHITE} ({nx} x {ny} x {nz}){RESET}")
    print(SEPARATOR)
    print()

# =====================================================
# Main
# =====================================================

def main():
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + " Lattice-Generator (CadQuery / STEP)" + RESET)
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
    cell_name, cell_step = chosen_variants["cell"]
    pad_name,  pad_step  = chosen_variants["pad"]
    hook_name, hook_step = chosen_variants["hook"]
    temp_name, temp_step = chosen_variants["temp"]
    bed_name,  bed_step  = chosen_variants["bed"]

    print()
    print(FG_WHITE + "Verwende Cell-Variante: " + FG_CYAN + (cell_name or "??") + RESET)
    print(FG_WHITE + "Cell-STEP: " + FG_CYAN + cell_step + RESET)

    print()
    print(FG_WHITE + "Pad-Variante:  " + FG_CYAN + (pad_name  or "keine") + RESET)
    print(FG_WHITE + "Hook-Variante: " + FG_CYAN + (hook_name or "keine") + RESET)
    print(FG_WHITE + "Temp-Variante: " + FG_CYAN + (temp_name or "keine") + RESET)
    print(FG_WHITE + "Bed-Variante:  " + FG_CYAN + (bed_name  or "keine") + RESET)

    if not os.path.exists(cell_step):
        raise FileNotFoundError(f"Eingabedatei '{cell_step}' nicht gefunden.")

    # ------------------------------------------------------------------
    # 2) Anzahl Zellen pro Richtung
    # ------------------------------------------------------------------
    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + " Anzahl der Zellen pro Richtung" + RESET)
    print(SEPARATOR)

    nx = ask_int(FG_WHITE + "Anzahl Zellen in X" + RESET, DEFAULT_N_CELLS)
    ny = ask_int(FG_WHITE + "Anzahl Zellen in Y" + RESET, DEFAULT_N_CELLS)
    nz = ask_int(FG_WHITE + "Anzahl Zellen in Z" + RESET, DEFAULT_N_CELLS)

    print_overview(nx, ny, nz)

    # ------------------------------------------------------------------
    # 3) Export-Ordner + Basisname
    # ------------------------------------------------------------------
    export_dir = os.path.join(os.getcwd(), EXPORT_DIR_NAME)
    os.makedirs(export_dir, exist_ok=True)

    variant_suffix = "_".join(variant_name_parts)
    base_name = variant_suffix
    final_output_stl = os.path.join(export_dir, base_name + ".stl")

    # ------------------------------------------------------------------
    # 4) Cell-STEP laden und normalisieren
    # ------------------------------------------------------------------
    print()
    print(FG_WHITE + f"Lade Cell-STEP: {FG_CYAN}{cell_step}{RESET}")
    cell_shape_raw = load_step_shape(cell_step)

    print(FG_WHITE + "Normalisiere Zelle auf Ursprung und ermittle Zellgröße..." + RESET)
    cell_shape, cell_size = normalize_shape_to_origin(cell_shape_raw)
    sx, sy, sz = cell_size

    print(
        FG_WHITE
        + f" Zellgröße (Bounding Box): "
        + FG_CYAN
        + f"sx={sx:.3f} mm, sy={sy:.3f} mm, sz={sz:.3f} mm"
        + RESET
    )

    # ------------------------------------------------------------------
    # 5) Strebendicke -> Schrittweiten
    # ------------------------------------------------------------------
    print()
    print(SEPARATOR)
    print(FG_MAGENTA + BOLD + " Strebendicke & Zell-Schrittweiten" + RESET)
    print(SEPARATOR)

    # reale Strebendicke in mm
    strut_mm_user = ask_float(
        FG_WHITE + "Strebendicke in mm (Überlappungssteuerung):" + RESET,
        DEFAULT_STRUT_MM,
        min_value=0.0,
    )

    if strut_mm_user <= 0.0:
        # keine Überlappung -> Schritt = Zellkantenlänge
        step_x, step_y, step_z = sx, sy, sz
        print(FG_WHITE + " Kante an Kante (keine Überlappung)." + RESET)
        print(
            FG_WHITE
            + " Kachel-Schrittweiten: "
            + FG_CYAN
            + f"step_x={step_x:.3f}, step_y={step_y:.3f}, step_z={step_z:.3f} mm"
            + RESET
        )
    else:
        # Schritt = Zellkante - Strebendicke (im gleichen mm-Maßstab)
        step_x = max(sx - strut_mm_user, sx * 0.1)
        step_y = max(sy - strut_mm_user, sy * 0.1)
        step_z = max(sz - strut_mm_user, sz * 0.1)

        if (
            step_x != sx - strut_mm_user
            or step_y != sy - strut_mm_user
            or step_z != sz - strut_mm_user
        ):
            print(
                FG_YELLOW
                + "Hinweis: Strebendicke war größer als 90% der Zellkante – "
                  "Schrittweite wurde zur Sicherheit begrenzt."
                + RESET
            )

        print(
            FG_WHITE
            + " Eingestellte Strebendicke (real): "
            + FG_CYAN
            + f"{strut_mm_user:.4f} mm"
            + RESET
        )
        print(
            FG_WHITE
            + " Effektive Schrittweiten (Zellkante - Strebendicke):"
            + RESET
        )
        print(
            FG_CYAN
            + f" step_x={step_x:.4f}, step_y={step_y:.4f}, step_z={step_z:.4f} mm"
            + RESET
        )

    step = (step_x, step_y, step_z)

    # ------------------------------------------------------------------
    # 6) Tiling in CadQuery
    # ------------------------------------------------------------------
    print()
    print(FG_WHITE + f"Erzeuge Gitter: {FG_CYAN}{nx} x {ny} x {nz}" + RESET)
    lattice_shape = build_lattice_cell_tiling(cell_shape, step, nx, ny, nz)

    # ------------------------------------------------------------------
    # 7) Optional: Union mit pad/hook/temp/bed (falls vorhanden)
    # ------------------------------------------------------------------
    result_shape = lattice_shape

    for label in ["pad", "hook", "temp", "bed"]:
        name, path = chosen_variants[label]
        if path is not None:
            result_shape = union_optional_variant(result_shape, path, label)

    # ------------------------------------------------------------------
    # 8) Export als STL
    # ------------------------------------------------------------------
    print()
    print(SEPARATOR)
    print(FG_WHITE + "Exportiere finales STL..." + RESET)
    print(FG_WHITE + f" Ausgabe-Pfad: {FG_CYAN}{final_output_stl}{RESET}")

    start_time = time.time()
    exporters.export(result_shape, final_output_stl, exportType="STL")
    elapsed = time.time() - start_time

    print(FG_GREEN + f"STL-Export fertig in {elapsed:.1f}s." + RESET)

    print(SEPARATOR)
    print(FG_GREEN + BOLD + " Alles fertig." + RESET)
    print(SEPARATOR)
    if os.path.exists(final_output_stl):
        print(FG_WHITE + " Ausgabe-STL: " + FG_CYAN + final_output_stl + RESET)
    else:
        print(FG_WHITE + " Finale STL konnte NICHT erzeugt werden." + RESET)
    print(SEPARATOR)


if __name__ == "__main__":
    main()
