# ============================================================
# USER SETTINGS (change only these)
# ============================================================
VIEW_MODE   = 'XZ'        # '3D', 'XY', 'XZ', 'YZ'
PROJECTION  = 'PARALLEL'  # 'PARALLEL' or 'PERSPECTIVE'
SHOW_TRIAD  = True        # True/False
SHOW_LEGEND = True        # True/False

# Keep current viewport view/orientation/zoom instead of forcing VIEW_MODE
KEEP_CURRENT_VIEW = True

# --- Legend formatting ---
LEGEND_NUMBER_FORMAT = 'SCIENTIFIC'   # 'SCIENTIFIC' or 'FIXED'
LEGEND_DECIMAL_PLACES = 2     # digits after decimal in legend (e.g., 3, 5, 8)
LEGEND_FONT_SIZE_PT   = 24    # legend font size in points (240 = 24 pt in X11 font string)
USE_VERDANA           = True  # try Verdana first; fallback to Helvetica if unavailable
# ============================================================

# --- Triad position ---
TRIAD_POS_X = 19              # triad position X (pixel-like viewport coords)
TRIAD_POS_Y = 17              # triad position Y (pixel-like viewport coords)
# ============================================================


# ============================================================
# USER SETTINGS (viewport)
# ============================================================
FIT_TO_WINDOW = True     # True = fit model/results to viewport window
ZOOM_FACTOR   = 1.0      # 1.0 = no extra zoom (use >1 to zoom in after fit)
# ============================================================


# ============================================================
# SIMPLE EXPORT BLOCK (paste at end)
# ============================================================
EXPORT_TIFFS    = True
OUT_DIR        = r'c:\temp\abaqus_tiff'
FILE_PREFIX    = 'result'
TIFF_SIZE       = (4096, 2304)     # (width, height)
FIT_EACH_IMAGE = True

# What to export (edit this list only)
EXPORT_LIST = [
    'S_Mises', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23',
    'LE11', 'LE22', 'LE33', 'LE12', 'LE13', 'LE23',
    'E11',  'E22',  'E33',  'E12',  'E13',  'E23',
    'U_Mag', 'U1', 'U2', 'U3'
]
# ============================================================


from abaqus import session
from abaqusConstants import *

vp = session.viewports[session.currentViewportName]

# White background
session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='White')

# --- Projection ---
vp.view.setProjection(projection=PARALLEL if PROJECTION.upper() == 'PARALLEL' else PERSPECTIVE)

# --- View orientation (explicit, not presets) ---
mode = VIEW_MODE.upper()

def set_plane_view(plane):
    if plane == 'XY':      # look along +Z, Y up
        vp.view.setViewpoint(viewVector=(0, 0, 1), cameraUpVector=(0, 1, 0))
    elif plane == 'XZ':    # look along +Y, Z up
        vp.view.setViewpoint(viewVector=(0, 1, 0), cameraUpVector=(0, 0, 1))
    elif plane == 'YZ':    # look along +X, Z up
        vp.view.setViewpoint(viewVector=(1, 0, 0), cameraUpVector=(0, 0, 1))
    elif plane == '3D':    # iso
        vp.view.setViewpoint(viewVector=(1, 1, 1), cameraUpVector=(0, 0, 1))
    else:
        raise ValueError("VIEW_MODE must be one of: '3D', 'XY', 'XZ', 'YZ'")
    vp.view.fitView()

if not KEEP_CURRENT_VIEW:
    set_plane_view(mode)

# --- Ensure contours exist (so a contour legend exists) ---
try:
    vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))
except:
    pass

# --- General annotation options ---
vp.viewportAnnotationOptions.setValues(
    triad=ON if SHOW_TRIAD else OFF,
    legend=ON if SHOW_LEGEND else OFF,
    compass=OFF,
    title=OFF,
    state=OFF,
    annotations=OFF
)

# --- Legend formatting ---
legend_fmt = LEGEND_NUMBER_FORMAT.upper()
if legend_fmt == 'SCIENTIFIC':
    _legend_number_format = SCIENTIFIC
elif legend_fmt == 'FIXED':
    _legend_number_format = FIXED
else:
    raise ValueError("LEGEND_NUMBER_FORMAT must be 'SCIENTIFIC' or 'FIXED'")

# --- Legend formatting ---
vp.viewportAnnotationOptions.setValues(
    legendBox=OFF,
    legendNumberFormat=_legend_number_format,
    legendDecimalPlaces=LEGEND_DECIMAL_PLACES
)

# --- Legend font ---
# NOTE: in '-*-...-*-*-240-*-*' the "240" is decipoints = 24.0 pt (x10).
legend_font_size = int(LEGEND_FONT_SIZE_PT * 10)

try:
    if USE_VERDANA:
        vp.viewportAnnotationOptions.setValues(
            legendFont='-*-verdana-bold-r-normal-*-*-%d-*-*-p-*-*-*' % legend_font_size
        )
    else:
        raise Exception("Skip Verdana")
except:
    vp.viewportAnnotationOptions.setValues(
        legendFont='-*-helvetica-bold-r-normal-*-*-%d-*-*-p-*-*-*' % legend_font_size
    )

# Move triad nearer the mesh (use user parameters)
vp.viewportAnnotationOptions.setValues(triadPosition=(TRIAD_POS_X, TRIAD_POS_Y))

# ------------------------------------------------------------
# Fit model/results to the viewport window (optional zoom after)
# ------------------------------------------------------------
try:
    if not KEEP_CURRENT_VIEW:
        if FIT_TO_WINDOW:
            vp.view.fitView()               # fits to current viewport window
        if ZOOM_FACTOR != 1.0:
            vp.view.zoom(ZOOM_FACTOR)       # >1 zoom in, <1 zoom out
except:
    pass


import os

def _ensure_dir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def _set_primary_from_tag(tag):
    """
    tag examples: S11, S_Mises, LE11, E12, U1, U_Mag
    returns True if variable exists and was set
    """
    try:
        if tag.startswith('S_') and tag == 'S_Mises':
            vp.odbDisplay.setPrimaryVariable('S', INTEGRATION_POINT, (INVARIANT, 'Mises'))
        elif tag.startswith('S'):
            vp.odbDisplay.setPrimaryVariable('S', INTEGRATION_POINT, (COMPONENT, tag))
        elif tag.startswith('LE'):
            vp.odbDisplay.setPrimaryVariable('LE', INTEGRATION_POINT, (COMPONENT, tag))
        elif tag.startswith('E'):
            vp.odbDisplay.setPrimaryVariable('E', INTEGRATION_POINT, (COMPONENT, tag))
        elif tag.startswith('U_') and tag == 'U_Mag':
            vp.odbDisplay.setPrimaryVariable('U', NODAL, (INVARIANT, 'Magnitude'))
        elif tag.startswith('U'):
            vp.odbDisplay.setPrimaryVariable('U', NODAL, (COMPONENT, tag))
        else:
            return False

        vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))
        return True
    except:
        return False

def _export_tiff(fullpath_no_ext):
    try:
        if FIT_EACH_IMAGE and not KEEP_CURRENT_VIEW:
            vp.view.fitView()
        session.printOptions.setValues(rendition=COLOR)
        session.tiffOptions.setValues(imageSize=TIFF_SIZE)
        session.printToFile(
            fileName=fullpath_no_ext,
            format=TIFF,
            canvasObjects=(vp,)
        )
        # --- auto-crop background (requires Pillow) ---
        try:
            from PIL import Image, ImageChops
            tiff_path = fullpath_no_ext + '.tiff'
            im = Image.open(tiff_path).convert('RGB')
            bg = Image.new('RGB', im.size, (255, 255, 255))  # white background
            diff = ImageChops.difference(im, bg)
            bbox = diff.getbbox()
            if bbox:
                im.crop(bbox).save(tiff_path)
        except:
            pass
        return True
    except:
        return False

if EXPORT_TIFFS:
    _ensure_dir(OUT_DIR)
    for tag in EXPORT_LIST:
        if _set_primary_from_tag(tag):
            fname = os.path.join(OUT_DIR, '%s_%s' % (FILE_PREFIX, tag))
            _export_tiff(fname)
# ============================================================



from pathlib import Path
from PIL import Image

INPUT_DIR = Path(r"C:\temp\abaqus_tiff")
OUTPUT_DIR = INPUT_DIR / "cropped"
EXTENSIONS = {".tif", ".tiff"}

TOL = 15


def is_bg(pixel, bg, tol):
    return (
        abs(int(pixel[0]) - int(bg[0])) <= tol and
        abs(int(pixel[1]) - int(bg[1])) <= tol and
        abs(int(pixel[2]) - int(bg[2])) <= tol
    )


def crop_background(img, tol=TOL):
    img = img.convert("RGB")
    px = img.load()
    w, h = img.size
    bg = px[0, 0]
    left = 0
    while left < w:
        found = False
        for y in range(h):
            if not is_bg(px[left, y], bg, tol):
                found = True
                break
        if found:
            break
        left += 1
    right = w - 1
    while right >= 0:
        found = False
        for y in range(h):
            if not is_bg(px[right, y], bg, tol):
                found = True
                break
        if found:
            break
        right -= 1
    top = 0
    while top < h:
        found = False
        for x in range(w):
            if not is_bg(px[x, top], bg, tol):
                found = True
                break
        if found:
            break
        top += 1
    bottom = h - 1
    while bottom >= 0:
        found = False
        for x in range(w):
            if not is_bg(px[x, bottom], bg, tol):
                found = True
                break
        if found:
            break
        bottom -= 1
    if left < right and top < bottom:
        return img.crop((left, top, right + 1, bottom + 1))
    return img


def main():
    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"Input folder not found: {INPUT_DIR}")
    OUTPUT_DIR.mkdir(exist_ok=True)
    files = [f for f in INPUT_DIR.iterdir() if f.suffix.lower() in EXTENSIONS]
    if not files:
        print("No TIFF files found.")
        return
    for f in files:
        try:
            img = Image.open(f)
            cropped = crop_background(img, tol=TOL)
            out_path = OUTPUT_DIR / f.name
            cropped.save(out_path)
            print(f"Done: {f.name}")
        except Exception as e:
            print(f"Failed: {f.name} -> {e}")


if __name__ == "__main__":
    main()




