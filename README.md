# auto-tilt-pair
This is a code for fiducial marker tracking under arbitrary conditions.

we clearly divide the workflow of fiducial marker correspondence into two stages: (i) initial transformation determination and (ii) local correspondence refinement. In the first stage, we model the transform estimation as a correspondence pair inquiry and verification problem. The local geometric constraint and invariant features are used to reduce the complexity of the problem. In the second stage, we encode the geometry distribution of the fiducial markers by a weighted Gaussian mixture model and introduce drift parameters to correct the effects of beam-induced motion and sample deformation. Comprehensive experiments on real-world datasets show the efficiency and effectiveness of the proposed algorithm. Especially, the proposed two-stage algorithm is able to produce an accurate tracking within an average of $\leqslant100$ ms per image, even for micrographs with hundreds of fiducial markers, which makes the real-time ET data processing possible.

This software accepts the file format defined by markerauto (proposed in the "A novel fully automatic scheme for fiducial marker-based alignment in electron tomography"). To compile the file, cmake and opencv is necessary. opencv 2.3 is suggested, because the bin of markerauto is precompiled by opencv 2.3. the cmake will output "autopairs" and "drawmatch" if the compilation is successful.

mkdir build
cd ./build
cmake ..
make -j 2

The preprocessed files about Nitrosop2 and Nitrosop3 are available at https://drive.google.com/drive/folders/1X-9m_dDaT1U3sjNx8SdL15xy2wMOwhRF?usp=sharing
To test the code, please run the exe as the "cmd" files in the subdirectory (the markerauto commend is also given within the "cmd" files, the users may restart all the thing from markerauto).

We are trying to fully integrate the auto-tilt-pair technique into markerauto 1.65, interested users could write to hanrenmin@gmail.com for free try.

