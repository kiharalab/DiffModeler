function calc_jl(stp, endp, pos, data, fsiv)
    dtotal = 0.0
    pos2 = zeros(3)

    for xp in stp[1]: endp[1]
        rx = convert(Float32, xp) - pos[1]
        rx = rx * rx
        for yp in stp[2]: endp[2]
            ry = convert(Float32, yp) - pos[2]
            ry = ry * ry
            for zp in stp[3]: endp[3]
                rz = convert(Float32, zp) - pos[3]
                rz = rz * rz
                d2 = rx + ry + rz
                v = data[xp, yp, zp] * exp(-1.5 * d2 * fsiv)
                dtotal += v
                pos2[1] += v * xp
                pos2[2] += v * yp
                pos2[3] += v * zp
            end
        end
    end

    return dtotal, pos2
end


function res_vec_jl(src_width, src_orig, src_dim, src_data, dest_width, dest_orig, dest_dim, dreso)
    fs = (dreso / src_width) * 0.5
    fs = fs * fs
    fsiv = 1.0 / fs
    fmaxd = (dreso / src_width) * 2.0
    dest_vec = zeros(dest_dim, dest_dim, dest_dim, 3)
    dest_data = zeros(dest_dim, dest_dim, dest_dim)
    dest_ss_data = zeros(dest_dim, dest_dim, dest_dim, 4)


    for x in 1:dest_dim
        for y in 1:dest_dim
            for z in 1:dest_dim
                xyz_arr = [x, y, z]
                pos = (xyz_arr * dest_width + dest_orig - src_orig) / src_width

                # convert pos to int
                pos = trunc.(Int32, pos)

                if pos[1] < 1 || pos[2] < 1 || pos[3] < 1 || pos[1] > src_dim[1] || pos[2] > src_dim[2] || pos[3] > src_dim[3]
                    continue
                end

                if src_data[pos[1], pos[2], pos[3]] == 0
                    continue
                end

                # Start Point and End Point
                stp = pos .- fmaxd
                endp = pos .+ fmaxd .+ 1

                # Clamp
                stp = clamp.(stp, 1, src_dim)
                endp = clamp.(endp, 1, src_dim)

                stp = trunc.(Int32, stp)
                endp = trunc.(Int32, endp)

                dtotal, pos2 = calc_jl(stp, endp, pos, src_data, fsiv)

                if dtotal == 0
                    continue
                end

                dest_data[x, y, z] = dtotal

                pos2 = pos2 / dtotal

                tmpcd = pos2 - pos

                dvec = sqrt(tmpcd[1]^2 + tmpcd[2]^2 + tmpcd[3]^2)

                if dvec == 0
                    dvec = 1.0
                end

                dest_vec[x, y, z, :] = tmpcd ./ dvec

            end
        end
    end
    return dest_data, dest_vec
end