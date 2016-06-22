# Sequencing error correction on High Throughtput Short Read (HTSR) data with spectral alignment on CUDA.

Advanced Algorithm and Parallel Programming project, Politecnico di Milano
Author: Brizzolari Cecilia (852399), Corneo Andrea (849793)
Academic Year: 2015

###Structure

This program consists in two part:
  * The first one is `spectrum_main.c` that contains the preprocessing part, filtering the `*.seq` file retrieving only the full reads. It then computes the spectrum and compresses it in a Bloom Filter structure.
  * The second one is `cuda_main.cu` that contains the two main kernel for voting and correcting the reads.

###Compilation

You could compile `spectrum_main.c` easily with GCC, remeber to include `spooky.c`, `city.c` and `hashes.c` contained in the subfolder `BloomFilter`.
To compile `cuda_main.cu` you should use NVCC, for simplicity it includes all the hashing function, so it doesn't require any include.

###Copyright

All the files contained in this repository, that are not linked esplicitly to other copyright, are released under GNU GPL v3 licenses. We refer to all the above file as "The Program".

All the files part of The Program will contain this declaration.

    Copyright (c) 2016, Brizzolari Cecilia, Corneo Andrea

    We refer as "The Program" to all the files and documentation present
    in this repository that reports this declaration.

    This file is part of The Program

    The Program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

As already said external file present in this repository have their copyright in the source code, please refer to it.
