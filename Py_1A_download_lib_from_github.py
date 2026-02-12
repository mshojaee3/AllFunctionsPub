# ------------------------- 
# Example call 
# -------------------------
#Py_1A_download_lib_from_github("pbc2d_version1.py")
#
#
RAW_BASE = "https://raw.githubusercontent.com/mshojaee3/AllFunctionsPub/main/"
def Py_1A_download_lib_from_github(py_name):
    # folder of this script (fallback to cwd)
    try:
        d = os.path.dirname(os.path.abspath(__file__))
    except Exception:
        d = os.getcwd()
    # save locally with spaces removed (import-friendly)
    local_name = py_name.replace(' ', '')
    local_path = os.path.join(d, local_name)
    if os.path.isfile(local_path):
        return local_path
    url = RAW_BASE + py_name.replace(' ', '%20')
    try:
        try:
            from urllib.request import urlopen  # Py3
        except ImportError:
            from urllib2 import urlopen         # Py2
        data = urlopen(url).read()
        with open(local_path, 'wb') as f:
            f.write(data)
        return local_path
    except Exception as e:
        print("Could not download:", url, "|", str(e))
        return None

