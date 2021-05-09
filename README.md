# quickLD v0.2.3: High-performance Computation of Linkage Disequilibrium on CPUs and GPUs


quickLD (qLD) is a tool to calculate Linkage disequilibrium (the non-random association between alleles at different loci), with highly efficient CPU and GPU kernels that utilize dense linear algebra (DLA) operations.

## Getting Started

quickLD is developed and tested in Linux distributions and should only work in Linux environments. The instructions bellow will get you a copy of the project up and running on your local machine for development and testing purposes. See References for further information in the paper behind this software.

### Clone the repository and get into the main directory

```
git clone git@github.com:pephco/quickLD.git
cd quickLD
```

### TL;WR
It is suggested that you (at least) briefly view the rest of the Readme. If you just want to test the program before investing more time into it, you can use the following commands to install and test-run quickLD in blis kernel mode:

```
./install_blis
./config install blis
./test_run full blis
```

If you want to use a simple heatmap viewer for the results of quickLD, we provide a Python implementation.
To install the required packages for the viewer, a requirements file is included in the main directory.
The simplest way to use it is the command:

```
pip3 install -r heatmap_viewer_reqs.txt
```

which installs the requirements to your global python packages.

### Basic Execution
The most basic execution, for diploid data, utilises these commands (change the appropriate fields to your values):

```
./bin/qLD-parseVCF -input 'PATH-TO-VCF' \
                   -output 'NEW-DIR' \
		   -chrom 'CHROMOSOME' \
                   -size 10

./bin/qLD-compute -input 'NEW-DIR' \
                  -output 'OUTPUT-DIR' \
                  -blis
```

### Prerequisites

* [make](https://www.gnu.org/software/make/)
* [BLIS](https://github.com/flame/blis)
* [OpenCL](https://github.com/OpenCL)

If you want to use the heatmap viewer:
* [Python3](https://www.python.org/download/releases/3.0/)

Most Systems probably already include make, if yours does not, install the ```make``` package using the package manager of your distro.

For BLIS we use v0.1.8 and a custom kernel, all of which are included in our project. It installs locally in the main directory and should not conflict with existing installations. To install BLIS, simply run:

```
./install_blis
```

There is no "one way" to install OpenCL, so you are on your own here. To provide some help, here is some info on the official pages of the GPU manufacturers:

* [NVIDIA](https://developer.nvidia.com/opencl)
* [AMD](https://www.amd.com/en/support/kb/release-notes/amdgpu-installation)

Alternatively, you can search your distribution's packages for a standalone (NVIDIA/AMD) openCL package, if it exists.

You will also need to install OpenCL headers for linking. In most repositories the package is named ```opencl-headers```:

* [Ubuntu](https://packages.ubuntu.com/search?keywords=opencl-headers)
* [Debian](https://packages.debian.org/search?keywords=opencl-headers)
* [Arch](https://www.archlinux.org/packages/extra/any/opencl-headers)
* [Fedora](https://fedora.pkgs.org/30/fedora-armhfp/opencl-headers-2.2-4.20190205git49f07d3.fc30.noarch.rpm.html)
* [Centos 7](https://centos.pkgs.org/7/epel-x86_64/opencl-headers-2.2-1.20180306gite986688.el7.noarch.rpm.html)

### Installing quickLD

To install quickLD, you only need to run the included configuration script and follow the prompt:

```
./config
```

For a brief explanation on how to use the configuration script with arguments instead of the prompt: 

```
./config help
```

Finally, to clean the binary files and the generated objects, run:

```
./config clean
```

## Running the tests

To test the installation, we included a version of the COVID-19 dataset, processed using [RaiSD](https://github.com/alachins/raisd).
The following sections break down the execution of the testing script, providing some more information about the separate programs in the project. 

As with all the provided utilities, 

```
./test_run help
```

briefly explains the arguments that test_run script supports. Running 

```
./test_run
```

without arguments, should produce similar to the following results:

### qLD-parse-VCF

The first segment parses the VCF file and breaks it into smaller chunks. It also creates a header file with information about the parsing, needed later in the pipeline. An example of the expected output (edited for clarity):

```
./test_run vcf

Parsing VCF file with qLD-parse-VCF (needed only once per dataset)
-----------------------------------------------------------------

Directory /home/$USER/$PATH-TO-MAIN-DIR/test_dataset/input/covid_19 created.

Part size: 1 MB.
Creating Multiple Files
Creating Header: $FILE

Total snips: 29817, Snip Size: 16622 bits, Snips per file: 504

 0% -- Creating File: $FILE_1.vcf.gz
 1% -- Renaming $FILE_1.vcf.gz
 1% -- Creating File: $FILE_2.vcf.gz
 3% -- Renaming $FILE_2.vcf.gz
 3% -- Creating File: $FILE_3.vcf.gz
        .
        .
        .
 98% -- Creating File: $FILE_59.vcf.gz
 99% -- Renaming $FILE_59.vcf.gz
100% -- Creating File: $FILE_60.vcf.gz
100% -- Renaming $FILE_60.vcf.gz
29817 Snips written
chrom 504 16622 29817 1 29817
get lines: 0.180347s, create Parts: 5.414211s, Total Time: 5.594558s 

Operation finished
```

### qLD-parse-2MDF

Next in the pipeline, *qLD-parse-2MDF* gets a directory with parsed VCF files (previously created) and creates MDF files. These files contain the same information as the VCF files, but compressed into 64-bit values. This reduces disk IO overhead in the computation step, and should be preferred as input. It also copies the header file along with the new files, since we still need its information for the computation step. An example of the expected output (edited for clarity):

```
./test_run mdf

Creating MDF files with qLD-parse-2MDF (needed only once per dataset)
---------------------------------------------------------------------

$FILE_1.mdf
$FILE_2.mdf
  .
  .
  .
$FILE_59.mdf
$FILE_60.mdf

Operation finished
```

### qLD-compute

The whole process up until now was modifying the input, making it compatible with the next and last tool, *qLD-compute*. To conduct pairwise LD calculations between different files or within the same file, you need to call this program. In the testing script we run a calculation using all the sites (lines in the VCF), including all the samples (columns in the VCF), in the full region of our input file. An example of the expected output (edited for clarity):

```
./test_run compute gpu

Running Computations with the selected modes:
---------------------------------------------

qLD_Correlator running with arguments:
	-input /home/$USER/$PATH-TO-MAIN-DIR/test_dataset/input/MDF_covid_19
	-output /home/$USER/$PATH-TO-MAIN-DIR/test_dataset/output/covid_19_report_gpu.txt
	-r2limit 0.2
	-ploidy haploid
	-mdf
	-gpu
GPU kernel active
MDF input active
Finished in: 0.040s
Output saved in: /home/$USER/$PATH-TO-MAIN-DIR/test_dataset/output

Operation finished
```
If everything was installed correctly, you should see the correct outputs for all the supported modes. If something fails, try to re-install the module that failed (BLIS, OpenCL).

## Processing your VCF files

After running the tests successfully, you are ready to run different VCF files. Since this software is in alpha state, not many VCF formats have been tested, apart from simplified versions of real datasets, processed by RaiSD. To run any VCF file, you can use the following logic (same as above):

1. Parse the VCF file using *qLD-parse-VCF* and the appropriate arguments.
2. (Optional but recommended) Create MDF files from the parsed VCFs (step 1) using *qLD-parse-2MDF*.
3. Use *qLD-compute* to calculate LD scores.

Every program has a help menu with all the arguments it uses, with a brief explanation of each argument. To invoke the help menu, run from the main directory (which will produce the result shown):

```
./bin/qLD-parse-VCF -help

```
VCF_parser manual
-----------------
	-input       inputFile
	-output      outputDir
	-size        parts_size
	-Wmin        snip index
	-Wmax        snip index
	-posWmin     snip pos
	-posWmax     snip pos
	-inputList   inputFile
	-chrom       chromosome
	-toSingleOutput

Description:
	-input     <STRING>  Specifies the name of the input alignment file.
	-output    <STRING>  Specifies the path of the output alignment files.
	-size      <INT>     Specifies the size of the memory footprint of the output
	                     alignment files in MB, if toSingleOutput is set,
	                     size is not needed. Supported file formats: VCF.

	-Wmin      <INT>     index of the minimum snip to be included, minimum 1 (default)
	-Wmax      <INT>     index of the maximum snip to be included,
	                     maximum total-Snips (default)
	-posWmin   <INT>     pos of the minimum snip to be included, must be valid
	-posWmax   <INT>     pos of the maximum snip to be included, must be valid
	-inputList <STRING>  input text file with the pos to keep
	-chrom     <STRING>  Specifies the chromosome to be extracted from the original VCF
	-toSingleOutput      Used to generate a new VCF that is part of the input file,
	                     -Wmin and -Wmax mandatory with this command
```
./bin/qLD-parse-2MDF -help

qLD-parse-2MDF manual
---------------------
	-input input_Directory
	-output output_Directory
	-sampleList input_File
	-ploidy correct_ploidy
	-seed seed_number
	-impute

Description:
	-input <STRING>		Specifies the path of the input alignment parsed files
	-output <STRING>	Specifies the path of the output alignment directory.
	-sampleList <STRING>	txt file with Format:
				"sample1
				 sample2
				 sample3
				   ...
				 sampleN"
				Specifies the name of the file that includes
				a list of valid samples from the input
				that will be selected for processing

	-ploidy <STRING>  Supported ploidy types:
	                       haploid:           single digit snip: ex. '0'
	                       phased_diploid:    double digit snip: ex. "0|0"
	                       unphased_diploid:  double digit snip:	ex. "0/0"
	-seed <INT>		Sets the random seed.
	-impute			Enables the imputation of missing data.
```
```
./bin/qLD-compute -help

qLD-compute manual
------------------
	-input        first_input_Folder
	-input2       second_input_Folder
	-output       output_File
	-ploidy       {haploid/phased_diploid/unphased_diploid}
	-posWmin1     snip_pos
	-posWmax1     snip_pos
	-posWmin2     snip_pos
	-posWmax2     snip_pos
	-sampleList   input_File
	-sampleList2  input_File
	-inputList    input_File
	-r2limit      value
	-threads      value
	-seed         value
	-sorted
	-blis
	-gpu

Description
	-input       <STRING>  Specifies the path of the first input alignment parsed files
	-input2      <STRING>  Specifies the path of the second input alignment parsed files
	                       Optional, uses Position Window 2 if set, else
	                       Position Window 1 is used
	-output      <STRING>  Specifies the name of the output alignment file.
	-ploidy      <STRING>  Supported ploidy types:
	                       haploid:           single digit snip: ex. '0'
	                       phased_diploid:    double digit snip: ex. "0|0"
	                       unphased_diploid:  double digit snip:	ex. "0/0"
	-sampleList  <STRING>  txt file with Format:
	                       "sample1
	                        sample2
	                        sample3
	                        ...
	                        sampleN"
	                       Specifies the name of the file that includes
	                       a list of valid samples from first input
	                       that will be selected for processing
	-sampleList2 <STRING>  Specifies the name of the file that includes
	                       a list of valid samples from second input
	                       that will be selected for processing
	                       Optional, uses seperate files for each input if set,
	                       else uses the file from -sampleList for both inputs

	Position Window 1 (mutually exclusive with input list)
	------------------------------------------------------
	-posWmin1    <INT>     pos of the minimum snip to be included to window 1
	-posWmax1    <INT>     pos of the maximum snip to be included to window 1

	Position Window 2 (mutually exclusive with input list)
	Optional, uses Input 2 if set, else Input 1 is used
	------------------------------------------------------
	-posWmin2    <INT>     pos of the minimum snip to be included to window 2
	-posWmax2    <INT>     pos of the maximum snip to be included to window 2

	Input list (mutually exclusive with inputs and position windows)
	----------------------------------------------------------------
	-inputList   <STRING>  csv file with Format:
	                       "input,posWmin1,posWmin2,input2,posWmin2,posWmax2"
	                       If single input is needed, duplicate the first 3 args.
	-sorted                (requires input list) Sorts Input List

	-competing             Enables competing task queue as an optimization
	-r2limit     <FLOAT>   the lowest r2 value to be included in the results (default 0.2)
	-threads     <INT>     Number of threads to run in parallel.
	                       Suggested to use physical core number at max.
	                       On your system this would be 4.
	-seed        <INT>     Sets the random seed.	
        -blis                  Use the blis framework for calculations
	                       (Your System is eligible for blis)
	-gpu                   use the gpu for calculations
	                       (there has to be a gpu in the system)
	-mdf                   use preprocessed files as input
	                       (needs mdf parsing your input before using this feature)
	                       Recommended, perfomance greatly enhanced
```

## Missing-data handling

quickLD provides two strategies to handle missing data. By default,
sites with missing data are discarded. Alternatively, the -impute
command-line flag can be used to deploy a probabilistic missing-data
imputation method on a per-site basis.

## Built With

* [BLIS](https://github.com/flame/blis) - The framework behind our CPU kernel
* [OpenCL](https://www.khronos.org/opencl/) - The framework behind our GPU kernel

## Authors

* **[C. Theodoris](https://github.com/StrayLamb2)**

## License

This project is licensed under the GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details

## Known Issues/Bugs

* Loading files in *qLD-compute* is not yet memory optimized and should hog good chunks of memory in big files/large inputLists.
* Consecutive runs with the same output name are not overwritting the old data, neither prompt the user, so use different output names or (re)move the old files    before re-using the name.
* Not supported data could produce results out of range (>1). In this case we save those values as '123.456' in the report for easy tracking/removal.
* Not supported data could produce empty reports (from trimming all the sites) without giving an error message.
* Invalid input directories can cause unexpected behavior.

## Version History

* 0.1.0:    Initial Release
* 0.2.0:    Properly placed competing queue as an optimization flag. The new default parallelization method is to assign a region of tasks 
            in the queue to each thread.
            Made config script compatible with older bash versions.
            Minor fixes.
* 0.2.1:    Using syrk method for single region calculations. 
            Several fixes.
* 0.2.2:    Added random seed selection in parse-2MDF and compute.
            Added -impute option in parse-2MDF to manipulate missing data.
            More fixes.
* 0.2.3:    Added -chrom flag to parser for VCF inputs with multiple chromosomes.
            Minor fixes.
