#=
Copyright (c) 2015 Miles Lubin and Yury Dvorkin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

# code to generate B, the weighted Lapacian:

#=
B_{ij} = { -\beta_ij, for each line connecting buses (i,j)
           sum_k \beta_kj, i = j, for incoming lines (k,j)
           0 otherwise
         }
where \beta_ij is the line susceptance

=#


function genB(buses, lines)

    n = length(buses)
    B = zeros(n,n)

    for l in lines
        i = l.tail
        j = l.head
        B[i,j] = -l.y
        B[j,i] = -l.y
        B[i,i] += l.y
        B[j,j] += l.y
    end

    return B
end

# generate \hat B by removing the row and column corresponding to the reference bus
function remove_col_and_row(B,refbus)
    @assert size(B,1) == size(B,2)
    n = size(B,1)
    return B[1:n .!= refbus, 1:n .!= refbus]
end
