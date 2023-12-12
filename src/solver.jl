# splits the grid into 2 subgrids, wich can be parallised. returns 2 sets of indices
function split_grid(sys::IsiSys)::NTuple{2, Vector{CartesianIndex}}
	if !iseven(sys.L)
		error("N has to be even in order to split the grid. N is $(sys.L)")
	end

	all_indices = CartesianIndices(sys.grid) |> collect
	inds1 = similar(all_indices, sys.L^2 ÷ 2)
	inds2 = similar(all_indices, sys.L^2 ÷ 2)

	for i in axes(all_indices, 1), j in axes(all_indices, 2)
		if !iseven(i+j)
			println("$(i+j) is odd, so $(all_indices[i]) is put into inds1[$(i÷2+1)]")
			inds1[(i+j)÷2 + 1] = all_indices[i,j]
		else
			inds2[(i+j)÷2] = all_indices[i,j]
		end
	end
	return (inds1, inds2)
end


