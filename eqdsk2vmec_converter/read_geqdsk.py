from freeqdsk import read_geqdsk as freeqdsk_read_geqdsk

def read_geqdsk(filename):
    """
    Reads a GEQDSK file using the freeqdsk library.
    
    Args:
        filename (str): The path to the GEQDSK file.
        
    Returns:
        dict: A dictionary containing the GEQDSK data.
    """
    print(f"Reading GEQDSK file: {filename} using freeqdsk")
    gdata = freeqdsk_read_geqdsk(filename)
    return gdata


