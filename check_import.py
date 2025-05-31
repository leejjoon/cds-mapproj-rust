import sys
print("Python sys.path:")
for p in sys.path:
    print(f"  - {p}")

print("\nPYTHONPATH environment variable:")
import os
print(os.environ.get("PYTHONPATH"))

print("\nAttempting to import mapproj_py...")
try:
    import mapproj_py
    print("Successfully imported mapproj_py")
    print(f"mapproj_py location: {mapproj_py.__file__}")
    print(f"mapproj_py contents: {dir(mapproj_py)}")
except ImportError as e:
    print(f"Failed to import mapproj_py: {e}")
    import traceback
    traceback.print_exc()
