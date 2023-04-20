# environment setup notes

## AnnData Version
After the release of `anndata v0.8.0`, the disk format of `.h5ad` was changed. This resulted in AnnDatas that were written using `anndata` versions 0.8.0 and greater to be unreadable with older versions. This will result in errors like the following:

```
KeyError: 'dict'

anndata._io.utils.AnnDataReadError: Above error raised while reading key '/layers' of type <class 'h5py._hl.group.Group'> from /.
```

If this error occurs when the rules try to read in the provided `.h5ad` files, try updating the `anndata` and `scanpy` versions using either of the commands:
```
conda update -c conda-forge anndata scanpy
```
```
conda install -c conda-forge anndata=0.9.1 scanpy=1.9.1
```
This workflow has been tested to work using `anndata=0.9.1` and `scanpy=1.9.3` in addition to the versions specied in the `environment.yml` file.
