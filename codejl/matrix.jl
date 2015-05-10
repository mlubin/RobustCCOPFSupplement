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
