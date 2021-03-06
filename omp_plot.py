# Copyright (c) 2016, Los Alamos National Security, LLC
# All rights reserved.
# Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.
# Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import matplotlib.pyplot as plt

Time = [3.00277 , 1.62415 , 1.17352 , 0.912296 , 0.836073 , 0.768095 , 0.748521 , 0.713546 , 0.707184 , 0.714413 , 0.706138 , 0.719206 , 0.725228 , 0.715083 , 0.710973 , 0.745081 , 0.705365 , 0.721511 , 0.709054 , 0.739904 , 0.728212 , 0.784384 , 0.759263 , 0.75875 , 0.743516 , 0.756969 , 0.752923 , 0.762673 , 0.747198 , 0.751036 , 0.75543 , 0.798351 , 0.751355 , 0.743142 , 0.75727 , 0.758155 , 0.768544 , 0.767764 , 0.769831 , 0.779022]
 

Procs = list(range(1, 41))

plt.plot(Procs, Time)
plt.title("Execution time of omp.c on 8192x8192 for 256 iterations")
plt.ylabel("Time(s)")
plt.xlabel("Number of Threads")
plt.savefig("omp_plot.png")
