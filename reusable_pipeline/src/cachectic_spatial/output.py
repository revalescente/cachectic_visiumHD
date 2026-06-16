def write_sample(sdata, block, sample, output_path):
    path = output_path.format(block=block, sample=sample)
    sdata.write(path)
