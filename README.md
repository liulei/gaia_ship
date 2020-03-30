# GAIA SHiP (盖亚方舟): Star cluster Hunting Pipeline
### Identifying star clusters in Gaia archive using friend of friend (FoF) method and isochrone fitting. 

## Introduction
There are 3 directories in this repo. 

- `data/`:
    - `K13.dat`: Kharchenko et al. (2013)
    - `CG18.txt` and `CG18b.txt`: these two catalogues together compose of the CG18 catalog.
    - `Bica.txt`: Bica et al. (2019)
    - `cal_all.txt`: full version of Tab. 1 in the paper.
    - `cat_new.txt`: full version of Tab. 3 in the paper.
    - `group/`: cluter members of 56 cluster groups. 
    
    **Note**: `cat_all.txt`, `cat_new.txt` and `group/sc_groupXXXX.txt` can be loaded with `pandas`, e.g.:
    
        df = pd.read_csv("cat_all.txt", delim_whitespace = True, header = 0)
    
    - `fof/npy`: star members of 2443 star cluster candidates, in `.npy` format, can be loaded with `arr = np.load()`. **Note! These files might be demaged if you download them individually from this repo via web browser.** I strongly recommend you down the whole repo via `git clone` if you want to use the `.npy` format.
    - `fof/csv`: star members of 2443 star cluster candidates, in `.csv` format, can be loaded with `df = pd.read_csv()`. Safe for individual downling with web browser.
    - `isochrone/0-10.dat`: `.dat` files that contain the isochrone tables downloaded from Padova group web interface. 
    - `isochrone/Z_python3.npy.tar.gz`: the `Z.npy` data file. Provided in compressed format. Please uncompress it before loading: `tar zxvf Z_python3.npy.tar.gz`.
    - `fit_iso/`: isochrone fitting results of 2443 cluster candidates. For Python 3:` fit = np.load(name_fit_iso, allow_pickle=True, encoding='latin1').flat[0]`. For Python 2: `fit = np.load(name_fit_iso).flat[0]`.
         
- `figure/`:
    - `4panel/`: 4-panels of 2443 star cluster candidates. 
    
- `src/`: the GAIA SHiP pipeline, see below.
    
**Note**: 

- **Very important! Most of the programs are developed with Python 2.7.13 (provided by conda 4.5.8). If you use Python 3, please change `print " " ` to `print(" ")` and pay attention to the difference between `/` (float division) and `//` (integer division). **

- If you make use of SHiP in your work, we require that you quote the pipeline link `https://github.com/liulei/gaia_ship` and reference the following paper:

  - `Liu, Lei & Pang, Xiaoying, "A catalog of newly identified star clusters in GAIA DR2", 2019, ApJS, 245, 32, arXiv:1910.12600`

- According to feedbacks from colleagues, **star members in `.npy` format might be demaged and not be loadable with `np.load()` if you download them individually from this repo via web browser**. We strongly recommend you download the whole repo via `git clone` if you want to use the `.npy` format member list. For convenience, I have prepared the member list in csv format. This guarantees the safe downloading via web browser with a small loss of precision.

- The whole pipeline (including the data and figure) is as large as 3 GB, which is actually not easy to download from GitHub. For conveniece, I have prepared the `src.tar.gz`, just in case you are only interested with the code. If you want to run isochrone fitting, an extra file `Z_python.py.tar.gz` in `data/isochrone/` is required.

- Current SHiP pipeline includes the data preparation, FoF, isochrone fitting and classification parts, so that you may construct the same catalog presented in the above paper. The data visualization part is not provided, since the programs are not well documented and the writings are messy. However they are still available upon request.

- The FoF clustering part and isochrone fitting part are the core of this pipeline. For those want to integrate them into your own pipeline, more detailed instructions are provided for the corresponding file below (`isochrone.py` and `procfof.py`). 

- Due to the file size limitation set by GitHub (< 100 MB), `Z.npy` (~ 129 MB) cannot be uploaded directly. I provide the compressed version 'Z_python3.npy.tar.gz'. You may download it and uncompress it with `tar zxvf Z_python3.npy.tar.gz`. I only provide `Z.npy` for `GAIA` system. To generate for your own photometry system, prepare `.dat` files in Padova webpage, use `load_dat()` and `save_npy()` in `isochrone.py`. . 

- You may use `npy2csv.py` to convert the npy format member list of every individual SC candidate to csv format which is readable by topcat.

Feel free to contact me (`liulei@shao.ac.cn`) if you have any problem.

## `procgaia.py`

**Input**:

- Raw GAIA DR2 data (`GaiaSource_*_*.csv.gz` file).
- `ali-gaia_dr2_source.txt`: listing all `.csv.gz` files.

**Output**:

- `gaia_segxxxx.npy`: 200 files, stored in array of structure `dtype_gaia`.

**Description**:

- Retrieve position, proper motion, magnitude, error information from raw data, and store them in 200 `gaia_seg` files. See `dtype_gaia` for details.

## `procsel.py`

**Input**:

- `gaia_segXXXX.npy`: 200 files.

**Output**:

- `selXXXX.npy`: 200 files, stored in array of structure `dtype_select`, which consists of two parts:
    - `idx`: star index in `gaia_segXXXX.npy`
    - `param`: 5 parameters: l, b, parallax, pmra, pmdec
    
**Description**:

- Retrieve all stars that satisfy the following criterion:
    - mag_G < 18
    - 0.2 < parallax < 7.0 (in mas)
    - |pmra| < 30 and |pmdec| < 30 (in mas/year)
    - |b| < 25 degree
- Positions have been evolved from epoch 2015.5 to epoch 2000.0 according to proper motion.

## `partition.py`

**Input**:

- `selXXXX.npy`: 200 files.

**Output**:

- `gaia_partition.npy`

**Description**:

- This serial program equal partitions the 3D parameter space (l, b, parallax) continuously, such that each partition contains roughly the same number of stars. 
- The typical scale of star cluster (`l_sc`) and the dispersion of parallax (`sigma_plx`) determine the minimum size of the partition in the corresponding dimension. See the program for the values.

## `procpartition.py`

**Input**:

- `gaia_partition.npy`: partition table, generated by `partition.py`.
- `selXXXX.npy`

**Output**:

- `stars_in_segxxxx.npy`: 200 files, each row is a list that contains star indexes for the corresponding partition. 

**Description**:

- Assign stars to the corresponding partition according to partition table in `gaia_partition.npy`. 


## `procfof.py`

This pipeline is originally designed for the data processing of large number of particles. Therefore every step takes parallalization into account, which makes it complicated. Now I realize that most of cases, users only cares about a single sky area, which greatly simplifies the usage. Below is the intruction for FoF clustering usage for a single sky area.

- The FoF method can only identify cluster from a smooth background. Therefore it might not work well if your sky area is too large. 
- The entry for FoF clustering is `runfof()`. You may find instructions on how to prepare data in example function `load_stars_single()`. 
- Some variables for controlling the FoF clustering:
    - `self.w`: the weight of every quantity. Adjust them according to the feature of your data.
    - `self.b_fof`: the linking length. Usually 0.2 times the mean partile seperation. 
    - `n_star_min`: a group is output with at least this number of stellar numbers. Change it according to your usage. 
- The output of the `fof()` function is `ginfos` and `l_keys1`. See `output` below for instructions. For single sky area FoF clustering, you may care about the star member list of each group in `l_keys1`.

**Input**:

- `stars_in_segXXXX.npy`
- `selXXXX.npy`

**Output**:

- `ginfos_pXXXX.npy`: 2D floating array. Each row is a list that gives the basic information of a star cluster identified in this partition. 
    - Meaning of each element in the list: cluster star number, l, b, r\_max, pmra, pmdec, r\_pm, parallax, r\_parallax.
- `keys_sel_pXXXX.npy`: each row is a list that contains all keys for that cluster. Each key is a 64 bit integer:
    - bit 0 - 31: star index in `selXXXX.npy`
    - bit 32 - 63: file index for `selXXXX.npy`

**Description**:

- For each partition, the program identifies a series of star clusters using FoF algorithm, and save the basic info of the cluster and keys to the specific stars in `ginfos_pXXXX.npy` and `keys_sel_pXXXX.npy`, respectively.
- Note! Star clusters given in this step are just intermediate results for each partition. They need to be further merged by `mergefof_key.py` as final product.

## `mergefof_key.py`

**Input**:

- `ginfos_pXXXX.npy`
- `keys_sel_pXXXX.npy`

**Output**:

- `ginfos_merge.npy`: basic info of each cluster after merge
- `keys_seg_merge.npy`: keys (file and star index) of each cluster for `gaia_segXXXX.npy`

**Description**:

- This program merges cluster from different partitions. At present the merge is just based on the intersection of keys (`keys_sel_pXXXX.npy`) from two clusters. See `is_merge()` function for more details. 

## `prockeyseg2sc.py`

**Input**:

- `keys_seg_pXXXX.npy`
- `gaia_segXXXX.npy`

**Output**:

- `fof_scXXXX.npy`

**Description**:

- This program retrieve stars from `gaia_segXXXX.npy` according to keys of each cluster recorded in `keys_seg_pXXXX.npy`. 
- Further analysis of star clusters will be totally based on `fof_scXXXX.npy`.

## `prociso.py`

**Input**:

- `fof_scXXXX.npy`
- `Z.npy`: isochrones of multiple Z and ages

**Output**:

- `fit_iso_scXXXX.npy`

**Description**:

- This program reads `g` and `b-r` info from `fof_scXXXX.npy`, fits isochrones to derive Z and age. The fitting is carried out by minizing $\bar{d^2} = \sum_{k = 1}^{n}(x_k - x_{k, nn})^2 / n$, where $x_k = (b-r, g)$ is color and magnitude of the $k~\mathrm{th}$ cluster member, $x_{k, nn}$ is the nearest neighbor of $k~\mathrm{th}$ member in the isochrone. 

## `isochrone.py`

This file provides all the necessary functions to for isochrone fitting. Below is a summary of its usage. 

- The main function for isochrone fitting is `fit_age_Z()`. An example is provided by fit():

        iso =   ISO()
        iso.load_npy('Z_python3.npy')
        d_fit   =   iso.fit_age_Z(g, b_r)
`g` and `b_r` are numpy arrays for photometry. `d_fit` is a dict that stores the fitting result. You may load it with `np.load(fit_iso_name, allow_pickle=True).flat[0]`. Pay attention to the extra `encoding='latin1'` param if you load our fitting results with Ptython 3.
 
- Originally I expect the users generate their `Z.npy` using `load_dat()` and `save_npy()` for their own photometry data. However it seems more convenienet to provide `Z.npy` directly. Please download it from `data/isochrone/Z_python3.npy.tar.gz` and uncompress it with `tar zxvf Z_python3.npy.tar.gz`. 

- As mentioned in the paper, there is a magnitude cut `g < 17` for bright stars. You may change or discard this criteria. 

- For efficiency, the fitting only deals with a maximum of 1000 stars, which is enough for most of star cluster cases. This is controlled by `nmax`. You may change or discard this criteria. 

- The fitting is only tested on apparent magnitude. You may try with absolute magnitude. If you succeed, please let me know. 

- You may want to plot the fitting result. I provide this function in `showisofit.py`. The input is `g`, `b_r` and the name of your fitting result file. Or just modify it to pass the fittign result dict. The problem is pretty simple. 

Below is old information on generating `Z.npy` for your own photometry system.

**Input**:

- `.dat` files that contain isochrone information downloaded from Padova group web interface (http://stev.oapd.inaf.it/cgi-bin/cmd). 

**Output**:

- `Z.npy`: metallicity and age table generated from `.dat` files.

**Description**:

The `ISO` class in this file provides the following functions:

- `load_dat()` and `save_npy()`: generate `Z.npy` from `.dat` files.


## `extract_sc_info.py`

**`Input`**:

- `fit_iso_scXXXX.npy`: isochrone fitting result ($t_\mathrm{age}$, $\bar{d^2}$)
- `fof_scXXXX.npy`: `g` and `b-r` for narrowness ($r_\mathrm{n}$) calculation, $n_{g<17}$

**`Output`**:

- `sc_info.txt`: id\_sc, ntot, $n_{g<17}$, $\bar{d^2}$, $r_\mathrm{n}$, $Z$ (logMsol), $t_\mathrm{age}$ (in Gyr), $\Delta_g$, $\Delta_{b-r}$, $d_\mathrm{plx}$, $d_\mathrm{iso}$, classification

**`Description`**:

This program reads sc info and classifies them.
  
## `match_K13.py`, `match_CG18.py` and `match_B19.py`

**`Input`**: 

- `ginfos_merge.npy`: position and radius of star clusters
- `sc_info.txt`: classification info
- `K13.dat` (K13, Kharchenko 2013)
- `CG18.dat` and `CG18b.dat` (CG18, Cantat-Gaudin 2018a,b)
- `B19.dat` (B19, Bica et al. 2019)

**`Output`**: cross matching catalog

    - col 1: line No. (0 indexed) in `ginfos_merge.npy`
    - col 2: position seperation in two catalogs
    - col 3: radius in `ginfos_merge.npy`
    - col 4: radius in reference 
    - col 5: line No. (0 indexed) in K13, CG18 and B19
    - col 6: classification

**`Description`**:

- This program cross matches the FoF sc of this work with K13,  CG18 and B19. Two clusters are regarded as matched if their distance is smaller than both of their radii. 

## `find_group.py`

**`Input`**:

- `cat_all.txt`

**`Output`**:

- `sc_groupXXXX.txt`: cluster members of a group.
- `sc_groups.txt`: number of clusters in each group. 

**Description**:

- This program reads all the class 1 cluster candidates, identifies cluster groups with standard FoF method, outputs groups that contain at least 2 members. 

**Note**:

- In this program and only in this program, we use the $d$, $l$, $b$ to $X$, $Y$, $Z$ conversion described by Eq. 3 of Conrad et al. (2017). This is different from the commonly used conversion adopted in our paper (Fig. 8). 
