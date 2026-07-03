import json
import math


def polygon_bounds(coordinates):
    ring = coordinates[0]
    xs = [point[0] for point in ring]
    ys = [point[1] for point in ring]
    return (
        math.floor(min(xs)),
        math.floor(min(ys)),
        math.ceil(max(xs)),
        math.ceil(max(ys)),
    )


def polygon_rotated_rectangle(coordinates):
    ring = coordinates[0]
    points = ring[:4]
    cx = sum(point[0] for point in points) / 4
    cy = sum(point[1] for point in points) / 4
    width = math.dist(points[0], points[1])
    height = math.dist(points[1], points[2])
    angle = math.degrees(math.atan2(points[1][1] - points[0][1], points[1][0] - points[0][0]))
    return {
        "center": [cx, cy],
        "width": width,
        "height": height,
        "angle_degrees": angle,
    }


def info_from_feature(feature, index=0):
    properties = feature.get("properties", {})
    coordinates = feature["geometry"]["coordinates"]
    x0, y0, x1, y1 = polygon_bounds(coordinates)
    sample = (
        properties.get("sample")
        or properties.get("name")
        or properties.get("sample_name")
        or f"sample_{index + 1}"
    )
    sample_key = properties.get("sample_key")
    if sample_key is None and properties.get("block") is not None:
        sample_key = f"{properties['block']}_{sample}"

    info = {
        "min_coordinate": [x0, y0],
        "max_coordinate": [x1, y1],
        "sample_key": sample_key or sample,
        "rotated_rectangle": polygon_rotated_rectangle(coordinates),
        "geometry": feature["geometry"],
        "properties": properties,
    }
    return sample, info


def load_samples(path, block):
    with open(path) as f:
        data = json.load(f)

    if data.get("type") == "FeatureCollection":
        samples = {}
        for index, feature in enumerate(data.get("features", [])):
            properties = feature.get("properties", {})
            if block is not None and properties.get("block") not in (None, block):
                continue
            sample, info = info_from_feature(feature, index=index)
            samples[sample] = info
        return samples

    return data[block]


def load_samples_by_block(path):
    with open(path) as f:
        data = json.load(f)

    if data.get("type") != "FeatureCollection":
        return data

    samples_by_block = {}
    for index, feature in enumerate(data.get("features", [])):
        properties = feature.get("properties", {})
        block = properties.get("block")
        if block is None:
            raise ValueError("GeoJSON sample feature is missing properties.block")
        sample, info = info_from_feature(feature, index=index)
        samples_by_block.setdefault(block, {})[sample] = info
    return samples_by_block


def rotated_rectangle(info):
    if "geometry" in info and info["geometry"].get("type") == "Polygon":
        return polygon_rotated_rectangle(info["geometry"]["coordinates"])

    if "rotated_rectangle" in info:
        return info["rotated_rectangle"]

    x0, y0 = info["min_coordinate"]
    x1, y1 = info["max_coordinate"]
    return {
        "center": [(x0 + x1) / 2, (y0 + y1) / 2],
        "width": x1 - x0,
        "height": y1 - y0,
        "angle_degrees": 0.0,
    }


def rotated_rectangle_corners(info):
    rect = rotated_rectangle(info)
    cx, cy = rect["center"]
    width = rect["width"]
    height = rect["height"]
    angle = math.radians(rect.get("angle_degrees", 0.0))
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)

    corners = []
    for dx, dy in (
        (-width / 2, -height / 2),
        (width / 2, -height / 2),
        (width / 2, height / 2),
        (-width / 2, height / 2),
    ):
        x = cx + dx * cos_a - dy * sin_a
        y = cy + dx * sin_a + dy * cos_a
        corners.append([x, y])
    return corners


def sample_bounds(info):
    if "geometry" in info and info["geometry"].get("type") == "Polygon":
        return polygon_bounds(info["geometry"]["coordinates"])

    if "rotated_rectangle" not in info:
        x0, y0 = info["min_coordinate"]
        x1, y1 = info["max_coordinate"]
        return int(x0), int(y0), int(x1), int(y1)

    corners = rotated_rectangle_corners(info)
    xs = [point[0] for point in corners]
    ys = [point[1] for point in corners]
    return (
        math.floor(min(xs)),
        math.floor(min(ys)),
        math.ceil(max(xs)),
        math.ceil(max(ys)),
    )
