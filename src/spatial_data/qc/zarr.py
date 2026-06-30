import matplotlib.pyplot as plt
from spatialdata import read_zarr


def guess_image_key(sdata):
    full = [key for key in sdata.images if "full" in key.lower()]
    if full:
        return full[0]
    hires = [key for key in sdata.images if "hires" in key.lower()]
    if hires:
        return hires[0]
    return next(iter(sdata.images))


def guess_labels_keys(sdata):
    keys = []
    for key in sdata.labels:
        lower = key.lower()
        if "nuclei" in lower or "fibres" in lower or "fiber" in lower:
            keys.append(key)
    return keys


def guess_shapes_keys(sdata):
    keys = []
    for key in sdata.shapes:
        lower = key.lower()
        if "roi" in lower or "interest" in lower:
            keys.append(key)
    return keys


def summarize_zarr(path):
    sdata = read_zarr(path)
    print("images:", list(sdata.images))
    print("labels:", list(sdata.labels))
    print("shapes:", list(sdata.shapes))
    print("tables:", list(sdata.tables))
    print("attrs:", dict(sdata.attrs))
    return sdata


def plot_zarr(
    path,
    coordinate_system=None,
    image_key=None,
    labels_keys=None,
    shapes_keys=None,
    output=None,
):
    import spatialdata_plot

    sdata = read_zarr(path)
    image_key = image_key or guess_image_key(sdata)
    labels_keys = guess_labels_keys(sdata) if labels_keys is None else labels_keys
    shapes_keys = guess_shapes_keys(sdata) if shapes_keys is None else shapes_keys

    plotter = getattr(sdata, "pl").render_images(image_key)
    for labels_key in labels_keys:
        plotter = plotter.pl.render_labels(labels_key, alpha=0.35)
    for shapes_key in shapes_keys:
        plotter = plotter.pl.render_shapes(
            shapes_key,
            fill_alpha=0.35,
            outline_width=1.0,
            outline_color="black",
            outline_alpha=1.0,
        )

    kwargs = {}
    if coordinate_system is not None:
        kwargs["coordinate_systems"] = coordinate_system
    plotter.pl.show(**kwargs)

    if output is not None:
        plt.savefig(output, dpi=160)

    return sdata
