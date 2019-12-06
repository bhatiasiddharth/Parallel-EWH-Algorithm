# Parallel EWH
This implementation is based on the following paper - [Parallelization of Error Weighted Hashing for Approximate k-Nearest neighbour search on GPU-CPU hybrid](https://www.comp.nus.edu.sg/~sbhatia/assets/pdf/bhatia2016parallelization.pdf). *Siddharth Bhatia, Mohan Pavan Kumar Badarla*.  IEEE International Conference on High Performance and Smart Computing, 2016.

Error Weighted Hashing (EWH) is a fast algorithm for Approximate k-Nearest neighbour search in Hamming space. It is more efficient than traditional Locality Sensitive Hashing algorithm (LSH) since it generates shorter list of strings for finding the exact distance from the query.

We have parallelized the EWH algorithm using GPU (Cuda) and using GPU-CPU hybrid (Cuda+OpenMP). 

## List of implementations
1. Sequential LSH
2. Distributed LSH (using MPI)
3. Sequential EWH
4. Parallel EWH using GPU (Cuda)
5. Parallel EWH using GPU-CPU hybrid (Cuda+OpenMP)

## Citation
If you use this code for your research, please consider citing our paper.

```
@inproceedings{bhatia2016parallelization,
  title={Parallelization of Error Weighted Hashing for Approximate k-Nearest neighbour search on GPU-CPU hybrid},
  author={Bhatia, Siddharth and Badarla, Mohan Pavan Kumar},
  booktitle={2016 IEEE 2nd International Conference on Big Data Security on Cloud (BigDataSecurity), IEEE International Conference on High Performance and Smart Computing (HPSC), and IEEE International Conference on Intelligent Data and Security (IDS)},
  pages={203--207},
  year={2016},
  organization={IEEE}
}
```
