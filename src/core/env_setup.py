"""
Environment bootstrapper for GROMACS GUI.
Handles headless Linux display/OpenGL dynamic library requirements transparently.
"""

import os
import sys
import subprocess
import ctypes


def prepare_headless_env():
    """Ensure PyQt6 loads without issues in headless Linux environments."""
    if sys.platform != "linux":
        return

    # Check if a display server is present
    if os.environ.get("DISPLAY"):
        return

    # Default QT QPA platform for headless environments
    if "QT_QPA_PLATFORM" not in os.environ:
        os.environ["QT_QPA_PLATFORM"] = "minimal"

    libs = [
        "libGL.so.1",
        "libxkbcommon.so.0",
        "libEGL.so.1",
        "libdbus-1.so.3",
        "libpcsclite.so.1",
        "libwayland-client.so.0",
        "libxcb-cursor.so.0",
    ]

    fake_dir = "/tmp/fake_libs"
    os.makedirs(fake_dir, exist_ok=True)

    missing = False
    for lib in libs:
        try:
            ctypes.CDLL(lib)
        except OSError:
            missing = True
            break

    if missing:
        # Build missing stubs if not present
        if not os.path.exists(os.path.join(fake_dir, "libGL.so.1")):
            dummy_c = os.path.join(fake_dir, "dummy.c")
            with open(dummy_c, "w") as f:
                f.write("void dummy_func() {}\n")

            for lib in libs:
                target = os.path.join(fake_dir, lib)
                if lib == "libdbus-1.so.3":
                    try:
                        nm_out = subprocess.check_output(
                            ["nm", "-D", "--undefined-only", "/usr/local/lib/python3.11/dist-packages/PyQt6/Qt6/lib/libQt6DBus.so.6"],
                            text=True,
                        )
                        funcs = {line.strip().split()[-1].split("@")[0] for line in nm_out.splitlines() if "dbus_" in line}
                    except Exception:
                        funcs = {"dbus_server_get_address", "dbus_bus_get"}
                    dbus_c = os.path.join(fake_dir, "dbus.c")
                    with open(dbus_c, "w") as f:
                        f.write("\n".join(f"void {fn}() {{}}" for fn in funcs))
                    map_file = os.path.join(fake_dir, "dbus.map")
                    with open(map_file, "w") as f:
                        f.write("LIBDBUS_1_3 { global: *; };\n")
                    subprocess.run(
                        ["gcc", "-shared", "-fPIC", "-Wl,-soname,libdbus-1.so.3", f"-Wl,--version-script={map_file}", dbus_c, "-o", target],
                        check=False,
                    )

                elif lib == "libxkbcommon.so.0":
                    try:
                        nm_out = subprocess.check_output(
                            ["nm", "-D", "--undefined-only", "/usr/local/lib/python3.11/dist-packages/PyQt6/Qt6/lib/libQt6Gui.so.6"],
                            text=True,
                        )
                        funcs = {line.strip().split()[-1].split("@")[0] for line in nm_out.splitlines() if "xkb_" in line}
                    except Exception:
                        funcs = {"xkb_context_new", "xkb_state_unref"}
                    xkb_c = os.path.join(fake_dir, "xkb.c")
                    with open(xkb_c, "w") as f:
                        f.write("\n".join(f"void {fn}() {{}}" for fn in funcs))
                    map_file = os.path.join(fake_dir, "xkb.map")
                    with open(map_file, "w") as f:
                        f.write("V_0.5.0 { global: *; };\n")
                    subprocess.run(
                        ["gcc", "-shared", "-fPIC", "-Wl,-soname,libxkbcommon.so.0", f"-Wl,--version-script={map_file}", xkb_c, "-o", target],
                        check=False,
                    )

                elif lib in ["libGL.so.1", "libEGL.so.1"]:
                    prefix = "gl" if "GL" in lib else "egl"
                    try:
                        nm_out = subprocess.check_output(
                            ["nm", "-D", "--undefined-only", "/usr/local/lib/python3.11/dist-packages/PyQt6/Qt6/lib/libQt6Gui.so.6"],
                            text=True,
                        )
                        funcs = set()
                        for line in nm_out.splitlines():
                            parts = line.strip().split()
                            if len(parts) >= 2 and parts[0] == "U":
                                s = parts[1].split("@")[0]
                                if s.lower().startswith(prefix):
                                    funcs.add(s)
                    except Exception:
                        funcs = {"glClear", "eglQueryString"}
                    c_file = os.path.join(fake_dir, f"{prefix}.c")
                    with open(c_file, "w") as f:
                        f.write("\n".join(f"void {fn}() {{}}" for fn in funcs))
                    subprocess.run(
                        ["gcc", "-shared", "-fPIC", f"-Wl,-soname,{lib}", c_file, "-o", target],
                        check=False,
                    )
                else:
                    subprocess.run(
                        ["gcc", "-shared", "-fPIC", f"-Wl,-soname,{lib}", dummy_c, "-o", target],
                        check=False,
                    )

        curr_ld = os.environ.get("LD_LIBRARY_PATH", "")
        if fake_dir not in curr_ld:
            os.environ["LD_LIBRARY_PATH"] = f"{fake_dir}:{curr_ld}" if curr_ld else fake_dir
            if "GROMACS_GUI_REEXEC" not in os.environ:
                os.environ["GROMACS_GUI_REEXEC"] = "1"
                os.execve(sys.executable, [sys.executable] + sys.argv, os.environ)
