import json

def read_from_json(filename):
    """Read from a given json file."""
    with open(filename) as json_file:
        data = json.load(json_file)

    return data
