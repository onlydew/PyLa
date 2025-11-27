import bpy
import sys
import os

def main():
    print("=== Blender Boolean Union Script ===")

    # -------------------------------
    # 1) Argumente nach "--" auslesen
    # -------------------------------
    argv = sys.argv
    print("sys.argv:", argv)

    if "--" not in argv:
        print("Fehler: Keine '--' in sys.argv gefunden. Erwartet: blender ... -P script.py -- lattice.stl pads.stl output.stl")
        sys.exit(1)

    idx = argv.index("--") + 1
    if len(argv) < idx + 3:
        print("Fehler: Zu wenige Argumente nach '--'.")
        print("Erwartet: lattice_stl pads_stl output_stl")
        sys.exit(1)

    lattice_path = os.path.abspath(argv[idx])
    pads_path    = os.path.abspath(argv[idx + 1])
    output_path  = os.path.abspath(argv[idx + 2])

    print("Lattice STL:", lattice_path)
    print("Pads STL   :", pads_path)
    print("Output STL :", output_path)

    if not os.path.exists(lattice_path):
        print("Fehler: Lattice-STL existiert nicht!")
        sys.exit(1)
    if not os.path.exists(pads_path):
        print("Fehler: Pad-STL existiert nicht!")
        sys.exit(1)

    # -------------------------------
    # 2) Szene leeren
    # -------------------------------
    print("Setze Blender auf Factory Settings zurück...")
    bpy.ops.wm.read_factory_settings(use_empty=True)

    # -------------------------------
    # 3) Lattice importieren (NEU: wm.stl_import)
    # -------------------------------
    print("Importiere Lattice-STL (wm.stl_import)...")
    res1 = bpy.ops.wm.stl_import(filepath=lattice_path)
    print("  Import Lattice Ergebnis:", res1)

    lattice_objs = list(bpy.context.selected_objects)
    if not lattice_objs:
        print("Fehler: Nach Lattice-Import wurden keine Objekte ausgewählt.")
        sys.exit(1)

    bpy.ops.object.select_all(action='DESELECT')
    for o in lattice_objs:
        o.select_set(True)
    bpy.context.view_layer.objects.active = lattice_objs[0]
    bpy.ops.object.join()
    lattice_obj = bpy.context.view_layer.objects.active
    lattice_obj.name = "LATTICE"

    print("Lattice-Objekt:", lattice_obj.name)

    # -------------------------------
    # 4) Pads importieren (NEU: wm.stl_import)
    # -------------------------------
    print("Importiere Pad-STL (wm.stl_import)...")
    bpy.ops.object.select_all(action='DESELECT')
    res2 = bpy.ops.wm.stl_import(filepath=pads_path)
    print("  Import Pads Ergebnis:", res2)

    pads_objs = list(bpy.context.selected_objects)
    if not pads_objs:
        print("Fehler: Nach Pad-Import wurden keine Objekte ausgewählt.")
        sys.exit(1)

    bpy.ops.object.select_all(action='DESELECT')
    for o in pads_objs:
        o.select_set(True)
    bpy.context.view_layer.objects.active = pads_objs[0]
    bpy.ops.object.join()
    pads_obj = bpy.context.view_layer.objects.active
    pads_obj.name = "PADS"

    print("Pads-Objekt:", pads_obj.name)

    # -------------------------------
    # 5) Boolean-Union anwenden
    # -------------------------------
    print("Erzeuge Boolean-Modifier (UNION)...")
    bpy.ops.object.select_all(action='DESELECT')
    lattice_obj.select_set(True)
    bpy.context.view_layer.objects.active = lattice_obj

    bool_mod = lattice_obj.modifiers.new(name="BoolUnion", type='BOOLEAN')
    bool_mod.operation = 'UNION'
    bool_mod.solver = 'EXACT'
    bool_mod.object = pads_obj

    print(" wende Boolean-Modifier an...")
    try:
        bpy.ops.object.modifier_apply(modifier=bool_mod.name)
    except Exception as e:
        print("Fehler beim Anwenden des Boolean-Modifiers:", e)
        sys.exit(1)

    # Pads löschen
    bpy.data.objects.remove(pads_obj, do_unlink=True)
    print("Pads-Objekt entfernt.")

    # -------------------------------
    # 6) Ergebnis exportieren (NEU: wm.stl_export)
    # -------------------------------
    print("Exportiere Ergebnis-STL (wm.stl_export)...")

    bpy.ops.object.select_all(action='DESELECT')
    lattice_obj.select_set(True)
    bpy.context.view_layer.objects.active = lattice_obj

    export_dir = os.path.dirname(output_path)
    if export_dir and not os.path.exists(export_dir):
        print("Erstelle Export-Verzeichnis:", export_dir)
        os.makedirs(export_dir, exist_ok=True)

    try:
        res_export = bpy.ops.wm.stl_export(
            filepath=output_path,
            ascii_format=False,
            export_selected_objects=True,   # nur das selektierte Bool-Ergebnis
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

    print("=== Blender Boolean Union Script fertig ===")

if __name__ == "__main__":
    main()
