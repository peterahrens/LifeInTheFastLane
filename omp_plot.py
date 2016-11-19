# Copyright (c) 2016, Peter Ahrens
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Los Alamos National Laboratories nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PETER AHRENS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# LA-CC-17-018

import matplotlib.pyplot as plt

Time = [3.00277 , 1.62415 , 1.17352 , 0.912296 , 0.836073 , 0.768095 , 0.748521 , 0.713546 , 0.707184 , 0.714413 , 0.706138 , 0.719206 , 0.725228 , 0.715083 , 0.710973 , 0.745081 , 0.705365 , 0.721511 , 0.709054 , 0.739904 , 0.728212 , 0.784384 , 0.759263 , 0.75875 , 0.743516 , 0.756969 , 0.752923 , 0.762673 , 0.747198 , 0.751036 , 0.75543 , 0.798351 , 0.751355 , 0.743142 , 0.75727 , 0.758155 , 0.768544 , 0.767764 , 0.769831 , 0.779022]
 

Procs = list(range(1, 41))

plt.plot(Procs, Time)
plt.title("Execution time of omp.c on 8192x8192 for 256 iterations")
plt.ylabel("Time(s)")
plt.xlabel("Number of Threads")
plt.savefig("omp_plot.png")
