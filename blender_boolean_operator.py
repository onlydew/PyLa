import bpy
import sys
import os


def import_and_join(path, label):
    """
    Importiert eine STL, joint alle importierten Objekte zu einem Mesh
    und gibt das resultierende Objekt zurück.
    """
    print(f"Importiere {label}-STL (wm.stl_import)...")
    res = bpy.ops.wm.stl_import(filepath=path)
    print(f"  Import {label} Ergebnis:", res)

    objs = list(bpy.context.selected_objects)
    if not objs:
        print(f"Fehler: Nach {label}-Import wurden keine Objekte ausgewählt.")
        sys.exit(1)

    bpy.ops.object.select_all(action='DESELECT')
    for o in objs:
        o.select_set(True)
    bpy.context.view_layer.objects.active = objs[0]
    bpy.ops.object.join()
    obj = bpy.context.view_layer.objects.active
    obj.name = label.upper()

    print(f"{label}-Objekt:", obj.name)
    return obj


def main():
    print("=== Blender Boolean Operator Script ===")

    argv = sys.argv
    print("sys.argv:", argv)

    if "--" not in argv:
        print("Fehler: Keine '--' in sys.argv gefunden.")
        print("Erwartet: blender ... -P script.py -- OP mesh_a.stl mesh_b.stl output.stl")
        sys.exit(1)

    idx = argv.index("--") + 1
    if len(argv) < idx + 4:
        print("Fehler: Zu wenige Argumente nach '--'.")
        print("Erwartet: OP mesh_a.stl mesh_b.stl output.stl")
        sys.exit(1)

    op_str     = argv[idx].upper()              # 'UNION' oder 'DIFFERENCE'
    mesh_a_path = os.path.abspath(argv[idx+1])
    mesh_b_path = os.path.abspath(argv[idx+2])
    output_path = os.path.abspath(argv[idx+3])

    if op_str not in {"UNION", "DIFFERENCE"}:
        print(f"Fehler: Ungültige Operation '{op_str}'. Erlaubt: UNION, DIFFERENCE.")
        sys.exit(1)

    print("Operation :", op_str)
    print("Mesh A STL:", mesh_a_path)
    print("Mesh B STL:", mesh_b_path)
    print("Output STL:", output_path)

    if not os.path.exists(mesh_a_path):
        print("Fehler: Mesh-A-STL existiert nicht!")
        sys.exit(1)
    if not os.path.exists(mesh_b_path):
        print("Fehler: Mesh-B-STL existiert nicht!")
        sys.exit(1)

    # -------------------------------
    # Szene leeren
    # -------------------------------
    print("Setze Blender auf Factory Settings zurück...")
    bpy.ops.wm.read_factory_settings(use_empty=True)

    # -------------------------------
    # Mesh A importieren
    # -------------------------------
    obj_a = import_and_join(mesh_a_path, "A")

    # -------------------------------
    # Mesh B importieren
    # -------------------------------
    obj_b = import_and_join(mesh_b_path, "B")

    # -------------------------------
    # Boolean anwenden (A (op) B)
    # -------------------------------
    print(f"Erzeuge Boolean-Modifier ({op_str})...")
    bpy.ops.object.select_all(action='DESELECT')
    obj_a.select_set(True)
    bpy.context.view_layer.objects.active = obj_a

    bool_mod = obj_a.modifiers.new(name="BoolOp", type='BOOLEAN')
    bool_mod.operation = op_str
    bool_mod.solver = 'EXACT'
    bool_mod.object = obj_b

    print("Wende Boolean-Modifier an...")
    try:
        bpy.ops.object.modifier_apply(modifier=bool_mod.name)
    except Exception as e:
        print("Fehler beim Anwenden des Boolean-Modifiers:", e)
        sys.exit(1)

    # Mesh B löschen
    bpy.data.objects.remove(obj_b, do_unlink=True)
    print("Mesh-B-Objekt entfernt.")

    # -------------------------------
    # Ergebnis exportieren
    # -------------------------------
    print("Exportiere Ergebnis-STL (wm.stl_export)...")

    bpy.ops.object.select_all(action='DESELECT')
    obj_a.select_set(True)
    bpy.context.view_layer.objects.active = obj_a

    export_dir = os.path.dirname(output_path)
    if export_dir and not os.path.exists(export_dir):
        print("Erstelle Export-Verzeichnis:", export_dir)
        os.makedirs(export_dir, exist_ok=True)

    try:
        res_export = bpy.ops.wm.stl_export(
            filepath=output_path,
            ascii_format=False,
            export_selected_objects=True,
            apply_modifiers=True,
            use_scene_unit=False,
        )
        print("  Export Ergebnis:", res_export)
    except Exception as e:
        print("Fehler beim Exportieren der STL:", e)
        sys.exit(1)

    if os.path.exists(output_path):
        print("Erfolgreich exportiert nach:", output_path)
    else:
        print("Export gemeldet, aber Datei existiert nicht?!")
        sys.exit(1)

    print("=== Blender Boolean Operator Script fertig ===")


if __name__ == "__main__":
    main()
