function ipSolve(grid::Array{Int,2}, optimizer)

    # Start a chronometer
    start = time()
    
    # Create the model
    model = Model(optimizer)

    m, n = size(grid)

    numberPos = Vector{Tuple{Int,Int}}()

    for i in 1:m
        for j in 1:n
            if grid[i,j] > 0
                push!(numberPos, (i, j))
            end
        end
    end

    p = length(numberPos)

    R = Vector{Vector{Array{Int,2}}}()
    for k = 1:p
        gridList = Vector{Array{Int,2}}()
        (i0, j0) = numberPos[k]
        val = grid[i0, j0]
        div = divisors(val)
        for (d1, d2) in vcat(div, reversePairs(div))
            if d1 > m || d2 > n
                continue
            end
            for i1 = 1:m - d1 + 1
                if !(i1 <= i0 < i1 + d1) # rectangle must contain its associated number
                    continue
                end
                for j1 = 1:n - d2 + 1
                    if !(j1 <= j0 < j1 + d2) # rectangle must contain its associated number
                        continue
                    end
                    # (i1, j1) are the coordinates of the top-left corner of the considered rectangle
                    rectGrid = Array{Int,2}(undef, m, n)
                    for i = 1:m
                        for j = 1:n
                            if i1 <= i < i1 + d1 && j1 <= j < j1 + d2
                                rectGrid[i,j] = 1
                            else
                                rectGrid[i,j] = 0
                            end
                        end
                    end
                    push!(gridList, rectGrid)
                end
            end
        end
        push!(R, gridList)
    end

    indices = variableIndices(R)

    @variable(model, x[indices], Bin)

    for k = 1:p
        (i0, j0) = numberPos[k]
        q = length(R[k])

        # Exactly one rectangle per number is active
        @constraint(model, sum(x[(k, l)] for l = 1:q) == 1)
    end

    # Rectangles do not overlap each other
    for i in 1:m,j in 1:n
        @constraint(model, sum(R[k][l][i, j] * x[(k, l)] for k = 1:p,l = 1:length(R[k])) == 1)
    end

    # @objective(model, Min, 0)

    # Solve the model
    optimize!(model)

    elapsedTime = time() - start

    status = JuMP.primal_status(model)
    if status != MOI.NO_SOLUTION
        sol = round.(Int, value.(x)) # solution of the LP problem

        # build a more compact representation of the solution
        rects = Vector{Tuple{Int,Int,Int,Int}}()
        for k = 1:p
            q = length(R[k])
            for l = 1:q
                if sol[(k, l)] == 1
                    i1, j1, i2, j2 = m + 1, n + 1, 0, 0
                    for i in 1:m, j in 1:n
                        if R[k][l][i,j] == 1
                            i1 = min(i1, i)
                            j1 = min(j1, j)
                            i2 = max(i2, i)
                            j2 = max(j2, j)
                        end
                    end
                    push!(rects, (i1, j1, i2, j2))
                    break
                end
            end
        end
        return true, elapsedTime, rects
    else
        return false, elapsedTime, Nothing
    end
end
