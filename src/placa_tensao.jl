function calcd2(pg, pf, n, nf, dad::Union{placa_fina}, pre)
    d2Rdx2,
    d2Rdxdy,
    d2Rdy2,
    d2Sdx2,
    d2Sdxdy,
    d2Sdy2,
    d3Rdx2dy,
    d3Rdx3,
    d3Rdxdy2,
    d3Rdy3,
    d3Sdx2dy,
    d3Sdx3,
    d3Sdxdy2,
    d3Sdy3,
    d4Rdx2dy2,
    d4Rdx3dy,
    d4Rdx4,
    d4Rdxdy3,
    d4Rdy4,
    d4Sdx2dy2,
    d4Sdx3dy,
    d4Sdx4,
    d4Sdxdy3,
    d4Sdy4,
    dRdx,
    dRdy,
    dSdx,
    dSdy,
    d5Sdx2dy3,
    d5Sdx3dy2,
    d5Sdx4dy,
    d5Sdx5,
    d5Sdxdy4,
    d5Sdy5,
    d3Rdx2dy,
    d3Rdx3,
    d3Rdxdy2,
    d3Rdy3,
    d4Rdx2dy2,
    d4Rdx3dy,
    d4Rdx4,
    d4Rdxdy3,
    d5Rdx2dy3,
    d5Rdx3dy2,
    d5Rdx4dy,
    d5Rdx5,
    d5Rdy5,
    d5Rdxdy4 = pre
    rs = pg - pf
    nx = n[1]
    ny = n[2]
    m1 = nf[1]
    m2 = nf[2]
    #dad.k.Distance of source and field points
    r1 = rs[1]
    r2 = rs[2]
    r = norm(rs)
    # Thin plate fundamental solutions
    theta = atan(r2, r1)

    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2

    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2

    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[1])

    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[2])

    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)

    a = 1


    for i = 1:2
        constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
        constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Rdx2[i] =
            2 * log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            )

        d2Rdxdy[i] =
            2 *
            dad.k.d[i] *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) -
            4 *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Rdy2[i] =
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) -
            8 *
            dad.k.d[i] *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        #--------------------------------------------------------------------------------
        # Cálculo das derivadas terceiras de R
        #--------------------------------------------------------------------------------

        # Cálculo da derivada terceira de R1 e R2 em relação a x conforme Shi e Bezine

        d3Rdx3[i] = 4 * (cos(theta) + dad.k.d[i] * sin(theta)) / (r * constA)

        # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
        d3Rdx2dy[i] =
            4 * (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) / (r * constA)

        # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
        d3Rdxdy2[i] =
            4 * (
                (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
            ) / (r * constA)

        # Cálculo da derivada terceira de R1 e R2 em relação a y conforme Shi e Bezine
        d3Rdy3[i] =
            4 * (
                dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
            ) / (r * constA)

        #--------------------------------------------------------------------------------
        # Cálculo das derivadas quartas de R
        #--------------------------------------------------------------------------------

        # Cálculo da derivada quarta de R1 e R2 em relação a x conforme Shi e Bezine
        d4Rdx4[i] =
            -4 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) /
            (r^2 * constA^2)

        # Cálculo da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
        d4Rdx3dy[i] =
            -4 / r^2 * (
                dad.k.d[i] / constA +
                2 * dad.k.e[i]^2 * sin(theta) * cos(theta) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )

        # Cálculo da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
        d4Rdx2dy2[i] =
            -4 / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 * dad.k.e[i]^2 * cos(theta)^2 / constA^2
            )

        # Cálculo da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
        d4Rdxdy3[i] =
            -4 / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 *
                dad.k.e[i]^2 *
                cos(theta) *
                (2 * dad.k.d[i] * cos(theta) + (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )

        # Cálculo da derivada quarta de R1 e R2 em relação a y conforme Shi e Bezine
        d4Rdy4[i] =
            -4 / r^2 * (
                (dad.k.d[i]^4 - dad.k.e[i]^4) / constA -
                2 *
                dad.k.e[i]^2 *
                cos(theta) *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                ) / constA^2
            )

        # Cálculo da derivada quinta de R1 e R2 em relação a x
        d5Rdx5[i] =
            (
                -4 *
                (cos(theta) + dad.k.d[i] * sin(theta)) *
                (
                    -1 - dad.k.d[i]^2 +
                    3 * dad.k.e[i]^2 +
                    (-1 + dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(2 * theta) -
                    2 * dad.k.d[i] * sin(2 * theta)
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de R1 e R2 em relação a x e y
        d5Rdx4dy[i] =
            (
                8 * (
                    dad.k.d[i] * cos(theta)^3 +
                    3 * (dad.k.d[i]^2 + dad.k.e[i]^2) * cos(theta)^2 * sin(theta) +
                    3 *
                    dad.k.d[i] *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    cos(theta) *
                    sin(theta)^2 +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de R1 e R2 em relação a x e y
        d5Rdx3dy2[i] =
            (
                8 * (
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta)^3 +
                    3 *
                    dad.k.d[i] *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    cos(theta)^2 *
                    sin(theta) +
                    3 * (dad.k.d[i]^2 + dad.k.e[i]^2)^2 * cos(theta) * sin(theta)^2 +
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2)^2 * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de R1 e R2 em relação a x e y
        d5Rdx2dy3[i] =
            (
                8 * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta)^3 +
                    3 * (dad.k.d[i]^4 - dad.k.e[i]^4) * cos(theta)^2 * sin(theta) +
                    3 *
                    dad.k.d[i] *
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^2 *
                    cos(theta) *
                    sin(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^3 * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de R1 e R2 em relação a x e y
        d5Rdxdy4[i] =
            (
                8 * (
                    (dad.k.d[i]^4 - 6 * dad.k.d[i]^2 * dad.k.e[i]^2 + dad.k.e[i]^4) *
                    cos(theta)^3 +
                    3 *
                    dad.k.d[i] *
                    (dad.k.d[i]^4 - 2 * dad.k.d[i]^2 * dad.k.e[i]^2 - 3 * dad.k.e[i]^4) *
                    cos(theta)^2 *
                    sin(theta) +
                    3 *
                    (dad.k.d[i]^2 - dad.k.e[i]^2) *
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^2 *
                    cos(theta) *
                    sin(theta)^2 +
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2)^3 * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de R1 e R2 em relação a y
        d5Rdy5[i] =
            (
                8 * (
                    dad.k.d[i] *
                    (dad.k.d[i]^4 - 10 * dad.k.d[i]^2 * dad.k.e[i]^2 + 5 * dad.k.e[i]^4) *
                    cos(theta)^3 +
                    3 *
                    (
                        dad.k.d[i]^6 - 5 * dad.k.d[i]^4 * dad.k.e[i]^2 -
                        5 * dad.k.d[i]^2 * dad.k.e[i]^4 + dad.k.e[i]^6
                    ) *
                    cos(theta)^2 *
                    sin(theta) +
                    3 *
                    dad.k.d[i] *
                    (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) *
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^2 *
                    cos(theta) *
                    sin(theta)^2 +
                    (dad.k.d[i]^2 - dad.k.e[i]^2) *
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^3 *
                    sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        #--------------------------------------------------------------------------------
        # Cálculo das derivadas segundas de S
        #--------------------------------------------------------------------------------
        d2Sdx2[i] =
            2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdxdy[i] =
            dad.k.e[i] * (log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            )) +
            2 *
            dad.k.d[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdy2[i] =
            2 *
            dad.k.d[i] *
            dad.k.e[i] *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) +
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        # Cálculo da derivada terceira de S1 e S2 em relação a x conforme Shi e Bezine
        d3Sdx3[i] = -2 * dad.k.e[i] * sin(theta) / (r * constA)

        # Cálculo da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
        d3Sdx2dy[i] = 2 * dad.k.e[i] * cos(theta) / (r * constA)

        # Cálculo da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
        d3Sdxdy2[i] =
            2 *
            dad.k.e[i] *
            (
                2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
            ) / (r * constA)

        # Cálculo da derivada terceira de S1 e S2 em relação a y conforme Shi e Bezine
        d3Sdy3[i] =
            2 *
            dad.k.e[i] *
            (
                (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
            ) / (r * constA)

        #--------------------------------------------------------------------------------
        # Cálculo das derivadas quartas de S
        #--------------------------------------------------------------------------------

        # Cálculo da derivada quarta de S1 e S2 em relação a x conforme Shi e Bezine
        d4Sdx4[i] =
            4 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) /
            (r^2 * constA^2)

        # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
        d4Sdx3dy[i] =
            2 * dad.k.e[i] / r^2 * (
                1 / constA -
                2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / constA^2
            )

        # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
        d4Sdx2dy2[i] =
            -4 *
            dad.k.e[i] *
            cos(theta) *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) / (r^2 * constA^2)

        # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
        d4Sdxdy3[i] =
            -2 * dad.k.e[i] / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                (
                    2 *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    cos(theta) *
                    (cos(theta) + dad.k.d[i] * sin(theta)) -
                    4 * dad.k.e[i]^2 * cos(theta)^2
                ) / constA^2
            )

        # Cálculo da derivada quarta de S1 e S2 em relação a y conforme Shi e Bezine
        d4Sdy4[i] =
            -4 * dad.k.e[i] / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                cos(theta) * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) / constA^2
            )

        #--------------------------------------------------------------------------------
        # Cálculo das derivadas quintas de S
        #--------------------------------------------------------------------------------

        # Cálculo da derivada quinta de S1 e S2 em relação a x
        d5Sdx5[i] =
            (
                4 *
                dad.k.e[i] *
                sin(theta) *
                (
                    -3 * cos(theta)^2 - 6 * dad.k.d[i] * cos(theta) * sin(theta) +
                    (-3 * dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de S1 e S2 em relação a x e y
        d5Sdx4dy[i] =
            (
                -4 *
                dad.k.e[i] *
                (
                    -cos(theta)^3 +
                    3 * (dad.k.d[i]^2 + dad.k.e[i]^2) * cos(theta) * sin(theta)^2 +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de S1 e S2 em relação a x e y
        d5Sdx3dy2[i] =
            (
                -4 *
                dad.k.e[i] *
                (
                    -2 * dad.k.d[i] * cos(theta)^3 -
                    3 * (dad.k.d[i]^2 + dad.k.e[i]^2) * cos(theta)^2 * sin(theta) +
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^2 * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de S1 e S2 em relação a x e y
        d5Sdx2dy3[i] =
            (
                4 *
                dad.k.e[i] *
                cos(theta) *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta)^2 +
                    3 *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    sin(theta) *
                    (
                        2 * dad.k.d[i] * cos(theta) +
                        (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    )
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de S1 e S2 em relação a x e y
        d5Sdxdy4[i] =
            (
                4 *
                dad.k.e[i] *
                (
                    4 * dad.k.d[i] * (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta)^3 +
                    3 *
                    (3 * dad.k.d[i]^4 + 2 * dad.k.d[i]^2 * dad.k.e[i]^2 - dad.k.e[i]^4) *
                    cos(theta)^2 *
                    sin(theta) +
                    6 *
                    dad.k.d[i] *
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^2 *
                    cos(theta) *
                    sin(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^3 * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

        # Cálculo da derivada quinta de S1 e S2 em relação a y
        d5Sdy5[i] =
            (
                4 *
                dad.k.e[i] *
                (
                    (5 * dad.k.d[i]^4 - 10 * dad.k.d[i]^2 * dad.k.e[i]^2 + dad.k.e[i]^4) *
                    cos(theta)^3 +
                    12 *
                    dad.k.d[i] *
                    (dad.k.d[i]^4 - dad.k.e[i]^4) *
                    cos(theta)^2 *
                    sin(theta) +
                    3 *
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) *
                    (dad.k.d[i]^2 + dad.k.e[i]^2)^2 *
                    cos(theta) *
                    sin(theta)^2 +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2)^3 * sin(theta)^3
                )
            ) / (
                r^3 *
                (
                    cos(theta)^2 +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)^2 +
                    dad.k.d[i] * sin(2 * theta)
                )^3
            )

    end

    #dad.k.Derivadas segundas da solução fundamental
    d2wdx2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
    d2wdxdy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
    d2wdy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))

    #dad.k.Derivadas terceiras da solução fundamental
    d3wdx3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
    d3wdx2dy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
    d3wdxdy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
    d3wdy3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))

    #dad.k.Derivadas quartas da solução fundamental
    d4wdx4 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
    d4wdx3dy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
    d4wdx2dy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
    d4wdxdy3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
    d4wdy4 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))

    #dad.k.Derivadas quintas da soluçao fundamental
    d5wdx5 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d5Rdx5[1] + C2 * d5Rdx5[2] + C3 * (d5Sdx5[1] - d5Sdx5[2]))
    d5wdx4dy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d5Rdx4dy[1] + C2 * d5Rdx4dy[2] + C3 * (d5Sdx4dy[1] - d5Sdx4dy[2]))
    d5wdx3dy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d5Rdx3dy2[1] + C2 * d5Rdx3dy2[2] + C3 * (d5Sdx3dy2[1] - d5Sdx3dy2[2]))
    d5wdx2dy3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d5Rdx2dy3[1] + C2 * d5Rdx2dy3[2] + C3 * (d5Sdx2dy3[1] - d5Sdx2dy3[2]))
    d5wdxdy4 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d5Rdxdy4[1] + C2 * d5Rdxdy4[2] + C3 * (d5Sdxdy4[1] - d5Sdxdy4[2]))
    d5wdy5 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d5Rdy5[1] + C2 * d5Rdy5[2] + C3 * (d5Sdy5[1] - d5Sdy5[2]))

    # Cálculo de f1; f2 e f3 conforme Shi e Bezine
    f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
    f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
    f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2

    # Cálculo de h1; h2; h3 e h4 conforme Shi e Bezine
    h1 = dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 - dad.k.D12 * nx * ny^2
    h2 =
        4 * dad.k.D16 * nx + dad.k.D12 * ny * (1 + nx^2) + 4 * dad.k.D66 * ny^3 -
        dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
    h3 =
        4 * dad.k.D26 * ny + dad.k.D12 * nx * (1 + ny^2) + 4 * dad.k.D66 * nx^3 -
        dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
    h4 = dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 - dad.k.D12 * nx^2 * ny


    #--------------Fernando------------------------------------
    # Cálculo das derivadas segundas de Mn
    d2mndx2 = -(f1 * d4wdx4 + f2 * d4wdx3dy + f3 * d4wdx2dy2)
    d2mndy2 = -(f1 * d4wdx2dy2 + f2 * d4wdxdy3 + f3 * d4wdy4)
    d2mndxdy = -(f1 * d4wdx3dy + f2 * d4wdx2dy2 + f3 * d4wdxdy3)

    # Cálculo das derivadas segundas de Vn

    d2vndx2 = -(h1 * d5wdx5 + h2 * d5wdx4dy + h3 * d5wdx3dy2 + h4 * d5wdx2dy3)
    d2vndy2 = -(h1 * d5wdx3dy2 + h2 * d5wdx2dy3 + h3 * d5wdxdy4 + h4 * d5wdy5)
    d2vndxdy = -(h1 * d5wdx4dy + h2 * d5wdx3dy2 + h3 * d5wdx2dy3 + h4 * d5wdxdy4)

    # Cálculo das segundas derivadas de Rci


    # Calculo das segundas derivadas de dwdn

    d3wdndx2 = (d3wdx3 * nx + d3wdx2dy * ny)
    d3wdndy2 = (d3wdxdy2 * nx + d3wdy3 * ny)
    d3wdndxdy = (d3wdx2dy * nx + d3wdxdy2 * ny)


    # Atribuição das variáveis da solução fundamental - atribuição de w_ij à variável
    # de saída w_est
    d2w = [
        d2wdx2 -d3wdndx2
        d2wdy2 -d3wdndy2
        d2wdxdy -d3wdndxdy
    ]


    # Atribuição das variáveis da derivada solução fundamental - atribuição de Vn_ij
    # à variável de saída Vn_est
    d2Vn = [
        d2vndx2 -d2mndx2
        d2vndy2 -d2mndy2
        d2vndxdy -d2mndxdy
    ]
    d2w, d2Vn
end

function calcdt(pg, pf, n, nf, dad::Union{placa_fina}, pre)
    R,
    S,
    d2Rdx2,
    d2Rdxdy,
    d2Rdy2,
    d2Sdx2,
    d2Sdxdy,
    d2Sdy2,
    d3Rdx2dy,
    d3Rdx3,
    d3Rdxdy2,
    d3Rdy3,
    d3Sdx2dy,
    d3Sdx3,
    d3Sdxdy2,
    d3Sdy3,
    d4Rdx2dy2,
    d4Rdx3dy,
    d4Rdx4,
    d4Rdxdy3,
    d4Rdy4,
    d4Sdx2dy2,
    d4Sdx3dy,
    d4Sdx4,
    d4Sdxdy3,
    d4Sdy4,
    dRdx,
    dRdy,
    dSdx,
    dSdy,
    extra = pre
    rs = pg - pf
    nx = n[1]
    ny = n[2]
    m1 = nf[1]
    m2 = nf[2]
    t1 = -m2
    t2 = m1
    #dad.k.Distance of source and field points
    r1 = rs[1]
    r2 = rs[2]
    r = norm(rs)
    # Thin plate fundamental solutions
    theta = atan(r2, r1)

    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2

    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2

    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[1])

    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[2])

    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)

    a = 1


    for i = 1:2
        constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
        constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        R[i] =
            r^2 *
            ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 3
            ) -
            4 *
            r^2 *
            dad.k.e[i] *
            sin(theta) *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        S[i] =
            r^2 *
            dad.k.e[i] *
            sin(theta) *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 3
            ) +
            r^2 *
            ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdx[i] =
            2 *
            r *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) -
            4 *
            r *
            dad.k.e[i] *
            sin(theta) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdy[i] =
            2 *
            r *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                dad.k.e[i]^2 * sin(theta)
            ) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) -
            4 *
            r *
            dad.k.e[i] *
            (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
            atan((dad.k.e[i] * sin(theta)), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdx[i] =
            r *
            dad.k.e[i] *
            sin(theta) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) +
            2 *
            r *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdy[i] =
            r *
            dad.k.e[i] *
            (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) +
            2 *
            r *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                dad.k.e[i]^2 * sin(theta)
            ) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))


        d2Rdx2[i] =
            2 * log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            )

        d2Rdxdy[i] =
            2 *
            dad.k.d[i] *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) -
            4 *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Rdy2[i] =
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) -
            8 *
            dad.k.d[i] *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdx2[i] =
            2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdxdy[i] =
            dad.k.e[i] * (log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            )) +
            2 *
            dad.k.d[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdy2[i] =
            2 *
            dad.k.d[i] *
            dad.k.e[i] *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) +
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d3Rdx3[i] = 4 * (cos(theta) + dad.k.d[i] * sin(theta)) / (r * constA)

        # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
        d3Rdx2dy[i] =
            4 * (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) / (r * constA)

        # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
        d3Rdxdy2[i] =
            4 * (
                (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
            ) / (r * constA)

        # Cálculo da derivada terceira de R1 e R2 em relação a y conforme Shi e Bezine
        d3Rdy3[i] =
            4 * (
                dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
            ) / (r * constA)
        # Cálculo da derivada terceira de S1 e S2 em relação a x conforme Shi e Bezine
        d3Sdx3[i] = -2 * dad.k.e[i] * sin(theta) / (r * constA)

        # Cálculo da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
        d3Sdx2dy[i] = 2 * dad.k.e[i] * cos(theta) / (r * constA)

        # Cálculo da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
        d3Sdxdy2[i] =
            2 *
            dad.k.e[i] *
            (
                2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
            ) / (r * constA)

        # Cálculo da derivada terceira de S1 e S2 em relação a y conforme Shi e Bezine
        d3Sdy3[i] =
            2 *
            dad.k.e[i] *
            (
                (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
            ) / (r * constA)

        #--------------------------------------------------------------------------------
        # Cálculo das derivadas quartas de S
        #--------------------------------------------------------------------------------

        # Cálculo da derivada quarta de S1 e S2 em relação a x conforme Shi e Bezine
        d4Sdx4[i] =
            4 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) /
            (r^2 * constA^2)

        # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
        d4Sdx3dy[i] =
            2 * dad.k.e[i] / r^2 * (
                1 / constA -
                2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / constA^2
            )

        # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
        d4Sdx2dy2[i] =
            -4 *
            dad.k.e[i] *
            cos(theta) *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) / (r^2 * constA^2)

        # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
        d4Sdxdy3[i] =
            -2 * dad.k.e[i] / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                (
                    2 *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    cos(theta) *
                    (cos(theta) + dad.k.d[i] * sin(theta)) -
                    4 * dad.k.e[i]^2 * cos(theta)^2
                ) / constA^2
            )

        # Cálculo da derivada quarta de S1 e S2 em relação a y conforme Shi e Bezine
        d4Sdy4[i] =
            -4 * dad.k.e[i] / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                cos(theta) * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) / constA^2
            )


    end
    dwdx =
        1 / (8 * pi * dad.k.D22) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
    dwdy =
        1 / (8 * pi * dad.k.D22) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))

    #dad.k.Derivadas segundas da solução fundamental
    d2wdx2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
    d2wdxdy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
    d2wdy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))

    #dad.k.Derivadas terceiras da solução fundamental
    d3wdx3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
    d3wdx2dy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
    d3wdxdy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
    d3wdy3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))

    #dad.k.Derivadas quartas da solução fundamental
    d4wdx4 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
    d4wdx3dy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
    d4wdx2dy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
    d4wdxdy3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
    d4wdy4 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))


    # Cálculo de f1; f2 e f3 conforme Shi e Bezine
    f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
    f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
    f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2

    # Cálculo de h1; h2; h3 e h4 conforme Shi e Bezine
    h1 = dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 - dad.k.D12 * nx * ny^2
    h2 =
        4 * dad.k.D16 * nx + dad.k.D12 * ny * (1 + nx^2) + 4 * dad.k.D66 * ny^3 -
        dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
    h3 =
        4 * dad.k.D26 * ny + dad.k.D12 * nx * (1 + ny^2) + 4 * dad.k.D66 * nx^3 -
        dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
    h4 = dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 - dad.k.D12 * nx^2 * ny

    dmndx = -(f1 * d3wdx3 + f2 * d3wdx2dy + f3 * d3wdxdy2)
    dmndy = -(f1 * d3wdx2dy + f2 * d3wdxdy2 + f3 * d3wdy3)
    dvndx = -(h1 * d4wdx4 + h2 * d4wdx3dy + h3 * d4wdx2dy2 + h4 * d4wdxdy3)
    dvndy = -(h1 * d4wdx3dy + h2 * d4wdx2dy2 + h3 * d4wdxdy3 + h4 * d4wdy4)

    dwdt = -(dwdx * t1 + dwdy * t2)
    d2wdndt = -(d2wdx2 * nx * t1 + d2wdxdy * (nx * t2 + ny * t1) + d2wdy2 * ny * t2)
    dmndt = -(dmndx * t1 + dmndy * t2)
    dvndt = -(dvndx * t1 + dvndy * t2)

    # @infiltrate
    [dwdt, -d2wdndt], [dvndt, -dmndt]
end

function int_dt_sing(x, dad, g_el, xi0)
    # Cálculo da distância do ponto fonte [xf,yf] ao ponto campo
    dx = (x[3, 1] - x[1, 1]) / 2 * 3
    dy = (x[3, 2] - x[1, 2]) / 2 * 3
    L = sqrt(dx^2 + dy^2)
    sx = dx / L
    sy = dy / L
    nx = sy
    ny = -sx
    # Ângulo entre r2 e r1
    theta = atan(dy, dx)
    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[1])
    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[2])
    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
    a = 1
    h_el = zeros(6)
    # g_el = zeros(6)
    # Integrate[N1d,(qsi,-1,1)]
    for x in [:NlogLqqo, :intNsr2, :intNsr, :intN]
        @eval $x = zeros(3)
    end
    intN[1] = 3 / 4
    intN[2] = 1 / 2
    intN[3] = intN[1]

    if xi0 ≈ 0
        #   2*Integrate[N1d[qsi] Log[L [Abs[qsi - qsio]]/2],(qsi,-1,1)]*Jac
        #  Jac=L/2
        NlogLqqo[1] = -(L * (1 + log(8) - 3 * log(L))) / 8
        NlogLqqo[2] = (L * (-3 + log(L / 2))) / 4
        NlogLqqo[3] = -(L * (1 + log(8) - 3 * log(L))) / 8
        intNsr[1] = -3 / 2
        intNsr[2] = 0
        intNsr[3] = 3 / 2
        intNsr2[1] = 9 / 4
        intNsr2[2] = -13 / 2
        intNsr2[3] = 9 / 4
    elseif abs(xi0) ≈ 2 / 3
        NlogLqqo[1] = (L * (-39 + 10 * log(5) - 27 * log(6) + 27 * log(L))) / 72
        NlogLqqo[2] = (L * (-5 + (25 * log(5)) / 6 - 3 * log(6) + 3 * log(L))) / 12
        NlogLqqo[3] = -(L * (3 + 25 * log(6 / 5) + log(36) - 27 * log(L))) / 72

        intNsr[1] = (3 * (-4 + (4 * log(5)) / 3)) / 4
        intNsr[2] = 3
        intNsr[3] = 0
        intNsr2[1] = (3 * (-9 / 5 - 3 * log(5))) / 4
        intNsr2[2] = (-9 + 6 * log(5)) / 2
        intNsr2[3] = (3 * (3 - log(5))) / 4
    else
        NlogLqqo[1] = (L * (-7 + 3 * log(L))) / 8
        NlogLqqo[2] = (L * log(L)) / 4
        NlogLqqo[3] = (L * (-1 + 3 * log(L))) / 8
        intNsr[1] = (3 * (-10 + 5 * log(2))) / 8
        intNsr[2] = 9 / 2 - (5 * log(8)) / 4
        intNsr[3] = (3 * (-2 + log(2))) / 8
        intNsr2[1] = 0
        intNsr2[2] = 0
        intNsr2[3] = 0
    end
    for x in [
        :d2Rdx2,
        :d2Rdx2ns,
        :d2Rdxdy,
        :d2Rdxdyns,
        :d2Rdxdys,
        :d2Rdy2,
        :d2Rdy2ns,
        :d2Rdy2s,
        :d2Sdx2,
        :d2Sdx2ns,
        :d2Sdx2s,
        :d2Sdxdy,
        :d2Sdxdyns,
        :d2Sdxdys,
        :d2Sdy2,
        :d2Sdy2ns,
        :d2Sdy2s,
        :d3Rdx2dy,
        :d3Rdx3,
        :d3Rdxdy2,
        :d3Rdy3,
        :d3Sdx2dy,
        :d3Sdx3,
        :d3Sdxdy2,
        :d3Sdy3,
        :d4Rdx2dy2,
        :d4Rdx3dy,
        :d4Rdx4,
        :d4Rdxdy3,
        :d4Rdy4,
        :d4Sdx2dy2,
        :d4Sdx3dy,
        :d4Sdx4,
        :d4Sdxdy3,
        :d4Sdy4,
        :d2Rdx2s,
    ]
        @eval $x = zeros(2)
    end
    for j = 1:3
        for i = 1:2
            constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
            constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            #----------------------------------------------------------------------
            # Cálculo das integrais das derivadas segundas de R
            #----------------------------------------------------------------------
            # Cálculo das integrais da derivada segunda de R1 e R2 em relação a x conforme Shi e
            # Bezine
            d2Rdx2s[i] = 4 * NlogLqqo[j]
            d2Rdx2ns[i] = 2 * log(constA / a^2) * intN[j] * L / 2
            d2Rdx2[i] = d2Rdx2s[i] + d2Rdx2ns[i]
            # Cálculo das integrais da derivada segunda de R1 e R2 em relação a x e y conforme Shi e Bezine
            d2Rdxdys[i] = 4 * dad.k.d[i] * NlogLqqo[j]
            d2Rdxdyns[i] =
                (2 * dad.k.d[i] * log(constA / a^2) - 4 * dad.k.e[i] * constB) *
                intN[j] *
                L / 2
            d2Rdxdy[i] = d2Rdxdys[i] + d2Rdxdyns[i]
            # Cálculo das integrais da derivada segunda de R1 e R2 em relação a y conforme Shi e Bezine
            d2Rdy2s[i] = 4 * (dad.k.d[i]^2 - dad.k.e[i]^2) * NlogLqqo[j]
            d2Rdy2ns[i] =
                (
                    2 * (dad.k.d[i]^2 - dad.k.e[i]^2) * log(a^2 * constA) -
                    8 * dad.k.d[i] * dad.k.e[i] * constB
                ) *
                intN[j] *
                L / 2
            d2Rdy2[i] = d2Rdy2s[i] + d2Rdy2ns[i]

            #--------------------------------------------------------------------------------
            # Cálculo das integrais das derivadas segundas de S
            #--------------------------------------------------------------------------------
            # Cálculo das integrais das derivadas segunda de S1 e S2 em relação a x conforme Shi e Bezine
            d2Sdx2s[i] = 0
            d2Sdx2ns[i] = 2 * constB * intN[j] * L / 2
            d2Sdx2[i] = d2Sdx2s[i] + d2Sdx2ns[i]
            # Cálculo das integrais da derivada segunda de S1 e S2 em relação a x e y conforme Shi e Bezine
            d2Sdxdys[i] = 2 * dad.k.e[i] * NlogLqqo[j]
            d2Sdxdyns[i] =
                (dad.k.e[i] * log(constA / a^2) + 2 * dad.k.d[i] * constB) * intN[j] * L / 2
            d2Sdxdy[i] = d2Sdxdys[i] + d2Sdxdyns[i]
            # Cálculo das integrais da derivada segunda de S1 e S2 em relação a y conforme Shi e Bezine
            d2Sdy2s[i] = 4 * dad.k.d[i] * dad.k.e[i] * NlogLqqo[j]
            d2Sdy2ns[i] =
                (
                    2 * dad.k.d[i] * dad.k.e[i] * log(constA / a^2) +
                    2 * (dad.k.d[i]^2 - dad.k.e[i]^2) * constB
                ) *
                intN[j] *
                L / 2
            d2Sdy2[i] = d2Sdy2s[i] + d2Sdy2ns[i]
            #--------------------------------------------------------------------------------
            # Cálculo das integrais das derivadas terceiras de R
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a x conforme Shi e Bezine
            d3Rdx3[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (cos(theta) + dad.k.d[i] * sin(theta)) /
                (constA)
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdx2dy[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (constA)
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdxdy2[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
                ) / (constA)
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a y conforme Shi e Bezine
            d3Rdy3[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) / (constA)
            #--------------------------------------------------------------------------------
            # Cálculo das derivadas quartas de R
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x conforme Shi e Bezine
            d4Rdx4[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) /
                (constA^2)
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdx3dy[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    dad.k.d[i] / constA +
                    2 * dad.k.e[i]^2 * sin(theta) * cos(theta) /
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    )^2
                )
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdx2dy2[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    (dad.k.d[i]^2 + dad.k.e[i]^2) / (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ) - 2 * dad.k.e[i]^2 * cos(theta)^2 / constA^2
                )
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdxdy3[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ) -
                    2 *
                    dad.k.e[i]^2 *
                    cos(theta) *
                    (
                        2 * dad.k.d[i] * cos(theta) +
                        (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) /
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    )^2
                )
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a y conforme Shi e Bezine
            d4Rdy4[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    (dad.k.d[i]^4 - dad.k.e[i]^4) / constA -
                    2 *
                    dad.k.e[i]^2 *
                    cos(theta) *
                    (
                        (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                        2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) / constA^2
                )

            #--------------------------------------------------------------------------------
            # Cálculo das das integrais derivadas terceiras de S
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a x conforme Shi e Bezine
            d3Sdx3[i] =
                -2 * (2 * intNsr[j] / L) * L / 2 * dad.k.e[i] * sin(theta) / (constA)
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
            d3Sdx2dy[i] =
                2 * (2 * intNsr[j] / L) * L / 2 * dad.k.e[i] * cos(theta) / (constA)
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
            d3Sdxdy2[i] =
                2 * (2 * intNsr[j] / L) * L / 2 *
                dad.k.e[i] *
                (
                    2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                ) / (constA)
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a y conforme Shi e Bezine
            d3Sdy3[i] =
                2 * (2 * intNsr[j] / L) * L / 2 *
                dad.k.e[i] *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                ) / (constA)
            #--------------------------------------------------------------------------------
            # Cálculo das derivadas quartas de S
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x conforme Shi e Bezine
            d4Sdx4[i] =
                4 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                sin(theta) *
                (cos(theta) + dad.k.d[i] * sin(theta)) / (constA^2)
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdx3dy[i] =
                2 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                (
                    1 / constA -
                    2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / constA^2
                )
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdx2dy2[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                cos(theta) *
                (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (constA^2)
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdxdy3[i] =
                -2 * dad.k.e[i] * (4 * intNsr2[j] / L^2) * L / 2 * (
                    (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                    (
                        2 *
                        (dad.k.d[i]^2 + dad.k.e[i]^2) *
                        cos(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        4 * dad.k.e[i]^2 * cos(theta)^2
                    ) / constA^2
                )
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a y conforme Shi e Bezine
            d4Sdy4[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                (
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                    cos(theta) * (
                        dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                        (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                    ) / constA^2
                )
        end
        # Integrais das derivadas segundas da solução fundamental
        d2wdx2 =
            1 / (8 * pi) * (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
        d2wdxdy =
            1 / (8 * pi) *
            (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
        d2wdy2 =
            1 / (8 * pi) * (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))
        # Integrais das derivadas terceira da solução fundamental
        d3wdx3 =
            1 / (8 * pi) * (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
        d3wdx2dy =
            1 / (8 * pi) *
            (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
        d3wdxdy2 =
            1 / (8 * pi) *
            (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
        d3wdy3 =
            1 / (8 * pi) * (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))
        # Integrais das derivadas quarta da solução fundamental
        d4wdx4 =
            1 / (8 * pi) * (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
        d4wdx3dy =
            1 / (8 * pi) *
            (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
        d4wdx2dy2 =
            1 / (8 * pi) *
            (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
        d4wdxdy3 =
            1 / (8 * pi) *
            (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
        d4wdy4 =
            1 / (8 * pi) * (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))
        # Cálculo de f1; f2 e f3 conforme Shi e Bezine
        f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
        f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
        f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2
        # Cálculo de h1; h2; h3 e h4 conforme Shi e Bezine
        h1 = dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 - dad.k.D12 * nx * ny^2
        h2 =
            4 * dad.k.D16 * nx + dad.k.D12 * ny * (1 + nx^2) + 4 * dad.k.D66 * ny^3 -
            dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
        h3 =
            4 * dad.k.D26 * ny + dad.k.D12 * nx * (1 + ny^2) + 4 * dad.k.D66 * nx^3 -
            dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
        h4 = dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 - dad.k.D12 * nx^2 * ny
        #--------------------------------------------------------------------------------
        # Cálculo das demais variáveis da solução fundamental
        #--------------------------------------------------------------------------------
        mn = -(f1 * d2wdx2 + f2 * d2wdxdy + f3 * d2wdy2)
        vn = -(h1 * d3wdx3 + h2 * d3wdx2dy + h3 * d3wdxdy2 + h4 * d3wdy3)
        # Cálculo das derivadas primeiras de Mn
        dmndx = -(f1 * d3wdx3 + f2 * d3wdx2dy + f3 * d3wdxdy2)
        dmndy = -(f1 * d3wdx2dy + f2 * d3wdxdy2 + f3 * d3wdy3)
        # Cálculo das derivadas primeiras de Vn
        dvndx = -(h1 * d4wdx4 + h2 * d4wdx3dy + h3 * d4wdx2dy2 + h4 * d4wdxdy3)
        dvndy = -(h1 * d4wdx3dy + h2 * d4wdx2dy2 + h3 * d4wdxdy3 + h4 * d4wdy4)
        #--------------------------------------------------------------------------------
        # Cálculo da derivada da solução fundamental
        #--------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------
        # Cálculo das demais variáveis da derivada da solução fundamental
        #--------------------------------------------------------------------------------

        t1 = sx
        t2 = sy


        d2wdndt =
            -(d2wdx2 * nx * t1 + d2wdxdy * (nx * t2 + ny * t1) + d2wdy2 * ny * t2) /
            dad.k.D22
        dmndt = -(dmndx * t1 + dmndy * t2) / dad.k.D22
        dvndt = -(dvndx * t1 + dvndy * t2) / dad.k.D22


        g32 = -d2wdndt

        h32 = -dmndt

        h31 = dvndt

        if xi0 ≈ -2 / 3
            h_el[2*j-1:2*j] = [h31 h32]
            g_el[2j] = g32
        elseif xi0 ≈ -1
            h_el[2*j-1:2*j] = [h31 h32]
        elseif xi0 ≈ 1
            h_el[7-2*j:8-2*j] = [h31 h32]
        else
            h_el[7-2*j:8-2*j] = [-h31 -h32]
            g_el[8-2j] = -g32
        end
        # @infiltrate
    end
    return g_el, h_el
end

function calc_dw_int(dad::placa_fina, forcas, desloc, npg = 8)
    nelem = size(dad.ELEM, 1)# Quantidade de elementos discretizados no contorno
    n_fis = size(dad.NOS, 1)
    n_internos = size(dad.pontos_internos, 1)
    n_cantos = size(dad.k.cantos, 1)
    H = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    G = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    q = zeros(2 * n_fis + n_cantos + n_internos)
    qsi, w = gausslegendre(npg)# Quadratura de gauss
    # qsi2, w2 = gausslegendre(2npg)# Quadratura de gauss
    # normal_fonte = calc_normais(dad)
    pre = [zeros(2) for idx = 1:50]

    d2w_int = zeros(n_internos, 3)

    for i = 1:n_internos

        pf = dad.pontos_internos[i, :] # Coordenada (x,y)
        d2w = zeros(3)
        nf = zeros(2)
        for j = 1:length(dad.ELEM)#Laço dos elementos
            elem_j = dad.ELEM[j]#Laço dos elementos
            x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
            Δelem = (x[end, :] - x[1, :]) # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet =
                (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                norm(Δelem)^2 + elem_j.ξs[1]
            N_geo, dN = calc_fforma(eet, elem_j)
            ps = N_geo' * x
            b = norm(ps' - pf) / norm(Δelem)
            eta, Jt = sinhtrans(qsi, eet, b)
            # eta,Jt=telles(qsi,eet)
            # Cálculo das integrais
            # intd2w, intd2Vn = integraelem_d2(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
            intd2w, intd2Vn = integraelem_d2(pf, nf, x, qsi, w, elem_j, dad, pre)

            int_gd2w = compute_q_d2(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
            inds = [2(elem_j.indices) .- 1 2(elem_j.indices)]'[:]
            # Calcula os deslocamentos [2x1] no nó interno
            d2w = d2w + intd2w * forcas[inds] - intd2Vn * desloc[inds] + int_gd2w
            # @infiltrate
        end

        d2Rc, d2wc = compute_d2Rw(pf, nf, dad, pre)

        d2w =
            d2w + d2wc * forcas[1+2n_fis+n_internos:end] -
            d2Rc * desloc[1+2n_fis+n_internos:end]
        d2w_int[i, :] = d2w
    end
    d2w_int
end
function calc_dwdt_cont(dad::placa_fina, forcas, desloc, npg = 8)
    nelem = size(dad.ELEM, 1)# Quantidade de elementos discretizados no contorno
    n_fis = size(dad.NOS, 1)
    n_internos = size(dad.pontos_internos, 1)
    n_cantos = size(dad.k.cantos, 1)
    H = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    G = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    q = zeros(2 * n_fis + n_cantos + n_internos)
    qsi, w = gausslegendre(npg)# Quadratura de gauss
    # qsi2, w2 = gausslegendre(2npg)# Quadratura de gauss
    normal_fonte = calc_normais(dad)
    pre = [zeros(2) for idx = 1:50]

    dwdt_cont = zeros(n_fis)

    for i = 1:n_fis
        d_c = 0

        pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
        nf = normal_fonte[i, :]
        caso = "contorno"

        for j = 1:length(dad.ELEM)#Laço dos elementos
            elem_j = dad.ELEM[j]#Laço dos elementos
            x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
            Δelem = (x[end, :] - x[1, :]) # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet =
                (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                norm(Δelem)^2 + elem_j.ξs[1]
            N_geo, dN = calc_fforma(eet, elem_j)
            ps = N_geo' * x
            b = norm(ps' - pf) / norm(Δelem)
            eta, Jt = sinhtrans(qsi, eet, b)
            # Cálculo das integrais
            intdwdt, intdqndt = integraelem_dt(pf, nf, x, qsi, w, elem_j, dad, pre)
            # @infiltrate
            qdt = compute_q_dt(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)

            if caso == "contorno"
                nosing = elem_j.indices .== i
                if sum(nosing) == 1
                    no_pf = findfirst(nosing)
                    xi0 = elem_j.ξs[no_pf]
                    # h, g = integraelemsing(pf, nf, x, qsi2, w2, elem_j, dad, xi0)
                    intdwdt, intdqndt = int_dt_sing(x, dad, intdwdt, xi0)
                end
            elseif caso == "canto"
                if j == dad.k.cantos[i-n_fis, 8]
                    xi0 = 1.0
                    intdwdt, intdqndt = int_dt_sing(x, dad, intdwdt, xi0)
                elseif j == dad.k.cantos[i-n_fis, 9]
                    xi0 = -1.0
                    intdwdt, intdqndt = int_dt_sing(x, dad, intdwdt, xi0)
                end
            end
            inds = [2(elem_j.indices) .- 1 2(elem_j.indices)]'[:]
            # Calcula os deslocamentos [2x1] no nó interno
            # @infiltrate
            # @show intdwdt
            d_c = d_c + intdwdt' * forcas[inds] - intdqndt' * desloc[inds] + qdt
        end

        Rt, wt = compute_dtRw(pf, nf, dad, pre)
        d_c =
            d_c + wt' * forcas[1+2n_fis+n_internos:end] -
            Rt' * desloc[1+2n_fis+n_internos:end]
        # @infiltrate
        dwdt_cont[i] = 2d_c
    end
    dwdt_cont
end
function integraelem_dt(pf, nf, x, eta, w, elem, dad::placa_fina, pre)
    intdwdt = zeros(Float64, 2 * size(elem))
    intdqndt = zeros(Float64, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        udt, pdt = calcdt(pg', pf, [sy, -sx], nf, dad, pre)
        # @infiltrate
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        # @infiltrate
        # @show k
        intdwdt += (udt' * Nm * dgamadqsi * w[k])'
        intdqndt += (pdt' * Nm * dgamadqsi * w[k])'

        # @infiltrate
    end
    # @show intdwdt
    intdwdt, intdqndt
end
function integraelem_d2(pf, nf, x, eta, w, elem, dad::placa_fina, pre)
    intd2w = zeros(Float64, 3, 2 * size(elem))
    intVn = zeros(Float64, 3, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        d2w, d2Vn = calcd2(pg', pf, [sy, -sx], nf, dad, pre)
        # @infiltrate
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        # @infiltrate
        intd2w += d2w * Nm * dgamadqsi * w[k]
        intVn += d2Vn * Nm * dgamadqsi * w[k]

        # @infiltrate
    end
    intd2w, intVn
end

function compute_q_d2(pf, nf, x, eta, Gauss_w, elem, dad, pre)

    # Inicialização da variável q_el
    q_el = zeros(3)
    # Início da integração numérica sobre o elemento
    for k = 1:size(Gauss_w, 1)# Integração da carga distribuída sobre o domínio
        N, dN = calc_fforma(eta[k], elem)
        pc = N' * x # Ponto de gauss interpolador
        rs = pc' - pf
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        n1 = sy
        n2 = -sx
        xf = pf[1]
        yf = pf[2]
        #dad.k.Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]
        r = norm(rs)
        rd1 = r1 / r
        rd2 = r2 / r
        m1 = nf[1]
        m2 = nf[2]
        theta = atan(r2, r1)
        nr = n1 * rd1 + n2 * rd2

        G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
        H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
        C1 =
            ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[1])
        C2 =
            ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[2])
        C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
        a = 1
        A = dad.k.carga[1]
        B = dad.k.carga[2]
        C = dad.k.carga[3]
        Clocal = A * xf + B * yf + C
        intd2Rdx2rdr,
        intd2Rdy2rdr,
        intd2Rdxdyrdr,
        intd2Rdx2ABrdr,
        intd2Rdy2ABrdr,
        intd2RdxdyABrdr,
        intd2Sdx2rdr,
        intd2Sdy2rdr,
        intd2Sdxdyrdr,
        intd2Sdx2ABrdr,
        intd2Sdy2ABrdr,
        intd2SdxdyABrdr,
        extra = pre

        for i = 1:2
            constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
            constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            intd2Rdx2rdr[i] = Clocal * (r^2 * (-1 + log((constA * r^2) / a^2)))
            intd2Rdy2rdr[i] =
                Clocal * (
                    r^2 * (
                        -dad.k.d[i]^2 - 4 * constB * dad.k.d[i] * dad.k.e[i] +
                        dad.k.e[i]^2 +
                        (dad.k.d[i]^2 - dad.k.e[i]^2) * log((constA * r^2) / a^2)
                    )
                )
            intd2Rdxdyrdr[i] =
                Clocal * (
                    r^2 * (
                        -dad.k.d[i] - 2 * constB * dad.k.e[i] +
                        dad.k.d[i] * log((constA * r^2) / a^2)
                    )
                )

            intd2Rdx2ABrdr[i] =
                (
                    2 *
                    r^3 *
                    (-2 + 3 * log((constA * r^2) / a^2)) *
                    (A * cos(theta) + B * sin(theta))
                ) / 9
            intd2Rdy2ABrdr[i] = 0
            intd2RdxdyABrdr[i] = 0

            intd2Sdx2rdr[i] = Clocal * (constB * r^2)
            intd2Sdy2rdr[i] =
                Clocal * (
                    r^2 * (
                        -(dad.k.d[i] * dad.k.e[i]) +
                        constB * (dad.k.d[i]^2 - dad.k.e[i]^2) +
                        dad.k.d[i] * dad.k.e[i] * log((constA * r^2) / a^2)
                    )
                )
            intd2Sdxdyrdr[i] =
                Clocal * (
                    (
                        r^2 * (
                            2 * constB * dad.k.d[i] - dad.k.e[i] +
                            dad.k.e[i] * log((constA * r^2) / a^2)
                        )
                    ) / 2
                )


            intd2Sdx2ABrdr[i] = 0
            intd2Sdy2ABrdr[i] = 0
            intd2SdxdyABrdr[i] = 0

        end
        int_d2Fnrrdx2 =
            1 / (8 * pi * dad.k.D22) * (
                C1 * (intd2Rdx2rdr[1] + intd2Rdx2ABrdr[1]) +
                C2 * (intd2Rdx2rdr[2] + intd2Rdx2ABrdr[2]) +
                C3 * (
                    (intd2Sdx2rdr[1] + intd2Sdx2ABrdr[1]) -
                    (intd2Sdx2rdr[2] + intd2Sdx2ABrdr[2])
                )
            )

        int_d2Fnrrdy2 =
            1 / (8 * pi * dad.k.D22) * (
                C1 * (intd2Rdy2rdr[1] + intd2Rdy2ABrdr[1]) +
                C2 * (intd2Rdy2rdr[2] + intd2Rdy2ABrdr[2]) +
                C3 * (
                    (intd2Sdy2rdr[1] + intd2Sdy2ABrdr[1]) -
                    (intd2Sdy2rdr[2] + intd2Sdy2ABrdr[2])
                )
            )

        int_d2Fnrrdxdy =
            1 / (8 * pi * dad.k.D22) * (
                C1 * (intd2Rdxdyrdr[1] + intd2RdxdyABrdr[1]) +
                C2 * (intd2Rdxdyrdr[2] + intd2RdxdyABrdr[2]) +
                C3 * (
                    (intd2Sdxdyrdr[1] + intd2SdxdyABrdr[1]) -
                    (intd2Sdxdyrdr[2] + intd2SdxdyABrdr[2])
                )
            )
        # Solução fundamental de deslocamento
        int_d2q = nr / r * [
            int_d2Fnrrdx2
            int_d2Fnrrdy2
            int_d2Fnrrdxdy
        ]
        # Cálculo dos componentes de q_el [2 X 2]
        q_el = q_el + int_d2q * dgamadqsi * Gauss_w[k]
    end
    #--------------------------------------------------------------------------------
    return q_el
end
function compute_q_dt(pf, nf, x, eta, Gauss_w, elem, dad, pre)

    # Inicialização da variável q_el
    q_el = 0
    # Início da integração numérica sobre o elemento
    for k = 1:size(Gauss_w, 1)# Integração da carga distribuída sobre o domínio
        N, dN = calc_fforma(eta[k], elem)
        pc = N' * x # Ponto de gauss interpolador
        rs = pc' - pf
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        n1 = sy
        n2 = -sx
        xf = pf[1]
        yf = pf[2]
        #dad.k.Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]
        r = norm(rs)
        rd1 = r1 / r
        rd2 = r2 / r
        m1 = nf[1]
        m2 = nf[2]
        t1 = -m2
        t2 = m1
        theta = atan(r2, r1)
        nr = n1 * rd1 + n2 * rd2

        G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
        H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
        C1 =
            ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[1])
        C2 =
            ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[2])
        C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
        a = 1
        A = dad.k.carga[1]
        B = dad.k.carga[2]
        C = dad.k.carga[3]
        Clocal = A * xf + B * yf + C
        intdRdxABrdr,
        intdRdxrdr,
        intdRdyABrdr,
        intdRdyrdr,
        intdSdxABrdr,
        intdSdxrdr,
        intdSdyABrdr,
        intdSdyrdr,
        intRABrdr,
        intRkrdr,
        intSABrdr,
        intSkrdr,
        extra = pre

        for i = 1:2


            intRkrdr[i] =
                Clocal * (
                    r^3 * (
                        -16 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        (
                            -7 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            -1 - dad.k.d[i]^2 +
                            dad.k.e[i]^2 +
                            (-1 + dad.k.d[i]^2 - dad.k.e[i]^2) * cos(2 * theta) -
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 16

            intSkrdr[i] =
                Clocal * (
                    r^3 * (
                        2 *
                        dad.k.e[i] *
                        (
                            -7 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) +
                        2 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            1 + dad.k.d[i]^2 - dad.k.e[i]^2 +
                            (1 - dad.k.d[i]^2 + dad.k.e[i]^2) * cos(2 * theta) +
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 16

            intRABrdr[i] =
                (
                    r^4 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        -40 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        (
                            -17 +
                            5 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            -1 - dad.k.d[i]^2 +
                            dad.k.e[i]^2 +
                            (-1 + dad.k.d[i]^2 - dad.k.e[i]^2) * cos(2 * theta) -
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 50

            intSABrdr[i] =
                (
                    r^4 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        2 *
                        dad.k.e[i] *
                        (
                            -17 +
                            5 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) +
                        5 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            1 + dad.k.d[i]^2 - dad.k.e[i]^2 +
                            (1 - dad.k.d[i]^2 + dad.k.e[i]^2) * cos(2 * theta) +
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 50

            intdRdxABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        -4 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) +
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 4

            intdRdyABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        -4 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + 2 * dad.k.d[i] * sin(theta)) +
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 4

            intdSdxABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        dad.k.e[i] *
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) +
                        4 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 8

            intdSdyABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        dad.k.e[i] *
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        (cos(theta) + 2 * dad.k.d[i] * sin(theta)) +
                        4 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 8

            intdRdxrdr[i] =
                Clocal * (
                    2 *
                    r^2 *
                    (
                        -6 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) +
                        (
                            -8 +
                            3 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 9

            intdRdyrdr[i] =
                Clocal * (
                    2 *
                    r^2 *
                    (
                        -6 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + 2 * dad.k.d[i] * sin(theta)) +
                        (
                            -8 +
                            3 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 9

            intdSdxrdr[i] =
                Clocal * (
                    r^2 * (
                        dad.k.e[i] *
                        (
                            -8 +
                            3 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) +
                        6 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 9

            intdSdyrdr[i] =
                -Clocal * (
                    r^2 * (
                        -(
                            dad.k.e[i] *
                            (
                                -8 +
                                3 * log(
                                    (
                                        r^2 * (
                                            dad.k.e[i]^2 * sin(theta)^2 +
                                            (cos(theta) + dad.k.d[i] * sin(theta))^2
                                        )
                                    ) / a^2,
                                )
                            ) *
                            (cos(theta) + 2 * dad.k.d[i] * sin(theta))
                        ) -
                        6 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 9

        end



        intdwdxrdr =
            1 / (8 * pi) * (
                C1 * (intdRdxrdr[1] + intdRdxABrdr[1]) +
                C2 * (intdRdxrdr[2] + intdRdxABrdr[2]) +
                C3 * ((intdSdxrdr[1] + intdSdxABrdr[1]) - (intdSdxrdr[2] + intdSdxABrdr[2]))
            )

        intdwdyrdr =
            1 / (8 * pi) * (
                C1 * (intdRdyrdr[1] + intdRdyABrdr[1]) +
                C2 * (intdRdyrdr[2] + intdRdyABrdr[2]) +
                C3 * ((intdSdyrdr[1] + intdSdyABrdr[1]) - (intdSdyrdr[2] + intdSdyABrdr[2]))
            )

        # int2 = -1*(intdwdxrdr*m1+intdwdyrdr*m2)*nr/dad.k.D22


        intdt = -1 * (intdwdxrdr * t1 + intdwdyrdr * t2) * nr / dad.k.D22
        # Cálculo dos componentes de q_el [2 X 2]
        # @infiltrate
        q_el = q_el + intdt * dgamadqsi * Gauss_w[k]
    end
    #--------------------------------------------------------------------------------
    return q_el
end
function compute_d2Rw(pf, nf, dad, pre)
    m1 = nf[1]
    m2 = nf[2]
    RS = zeros(3, size(dad.k.cantos, 1))
    WS = zeros(3, size(dad.k.cantos, 1))
    for i = 1:size(dad.k.cantos, 1)
        pc = dad.k.cantos[i, 2:3]
        na1 = dad.k.cantos[i, 4]
        na2 = dad.k.cantos[i, 5]
        nd1 = dad.k.cantos[i, 6]
        nd2 = dad.k.cantos[i, 7]
        rs = pc - pf
        #dad.k.Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]

        G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
        H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
        C1 =
            ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[1])
        C2 =
            ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[2])
        C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
        # a = constFS
        a = 1
        theta = atan(r2, r1)
        r = norm(rs)
        d2Rdxdy,
        d2Rdy2,
        d2Sdx2,
        d2Sdxdy,
        d2Sdy2,
        d4Rdx2dy2,
        d4Rdx3dy,
        d4Rdx4,
        d4Rdxdy3,
        d4Rdy4,
        d4Sdx2dy2,
        d4Sdx3dy,
        d4Sdx4,
        d4Sdxdy3,
        d4Sdy4,
        d2Rdx2,
        extra = pre
        for i = 1:2
            constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
            constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            d2Rdx2[i] =
                2 * log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                )

            ####

            d2Rdxdy[i] =
                2 *
                dad.k.d[i] *
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) -
                4 *
                dad.k.e[i] *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
            ###
            d2Rdy2[i] =
                2 *
                (dad.k.d[i]^2 - dad.k.e[i]^2) *
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) -
                8 *
                dad.k.d[i] *
                dad.k.e[i] *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            d3Rdx3[i] = 4 * (cos(theta) + dad.k.d[i] * sin(theta)) / (r * constA)

            # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdx2dy[i] =
                4 * (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (r * constA)

            # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdxdy2[i] =
                4 * (
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
                ) / (r * constA)

            # Cálculo da derivada terceira de R1 e R2 em relação a y conforme Shi e Bezine
            d3Rdy3[i] =
                4 * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) / (r * constA)
            # Cálculo da derivada quarta de R1 e R2 em relação a x conforme Shi e Bezine
            d4Rdx4[i] =
                -4 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) /
                (r^2 * constA^2)

            # Cálculo da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdx3dy[i] =
                -4 / r^2 * (
                    dad.k.d[i] / constA +
                    2 * dad.k.e[i]^2 * sin(theta) * cos(theta) /
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    )^2
                )

            # Cálculo da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdx2dy2[i] =
                -4 / r^2 * (
                    (dad.k.d[i]^2 + dad.k.e[i]^2) / (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ) - 2 * dad.k.e[i]^2 * cos(theta)^2 / constA^2
                )

            # Cálculo da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdxdy3[i] =
                -4 / r^2 * (
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ) -
                    2 *
                    dad.k.e[i]^2 *
                    cos(theta) *
                    (
                        2 * dad.k.d[i] * cos(theta) +
                        (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) /
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    )^2
                )

            # Cálculo da derivada quarta de R1 e R2 em relação a y conforme Shi e Bezine
            d4Rdy4[i] =
                -4 / r^2 * (
                    (dad.k.d[i]^4 - dad.k.e[i]^4) / constA -
                    2 *
                    dad.k.e[i]^2 *
                    cos(theta) *
                    (
                        (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                        2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) / constA^2
                )




            d2Sdx2[i] =
                2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            d2Sdxdy[i] =
                dad.k.e[i] * (log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                )) +
                2 *
                dad.k.d[i] *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))


            d2Sdy2[i] =
                2 *
                dad.k.d[i] *
                dad.k.e[i] *
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) +
                2 *
                (dad.k.d[i]^2 - dad.k.e[i]^2) *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            # Cálculo da derivada quarta de S1 e S2 em relação a x conforme Shi e Bezine
            d4Sdx4[i] =
                4 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) /
                (r^2 * constA^2)

            # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdx3dy[i] =
                2 * dad.k.e[i] / r^2 * (
                    1 / constA -
                    2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / constA^2
                )

            # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdx2dy2[i] =
                -4 *
                dad.k.e[i] *
                cos(theta) *
                (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (r^2 * constA^2)

            # Cálculo da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdxdy3[i] =
                -2 * dad.k.e[i] / r^2 * (
                    (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                    (
                        2 *
                        (dad.k.d[i]^2 + dad.k.e[i]^2) *
                        cos(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        4 * dad.k.e[i]^2 * cos(theta)^2
                    ) / constA^2
                )

            # Cálculo da derivada quarta de S1 e S2 em relação a y conforme Shi e Bezine
            d4Sdy4[i] =
                -4 * dad.k.e[i] / r^2 * (
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                    cos(theta) * (
                        dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                        (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                    ) / constA^2
                )
        end

        d2wdx2 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
        d2wdxdy =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
        d2wdy2 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))

        #dad.k.Derivadas quartas da solução fundamental
        d4wdx4 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
        d4wdx3dy =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
        d4wdx2dy2 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
        d4wdxdy3 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
        d4wdy4 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))


        g1a = (dad.k.D12 - dad.k.D11) * na1 * na2 + dad.k.D16 * (na1^2 - na2^2)
        g2a = 2 * (dad.k.D26 - dad.k.D16) * na1 * na2 + 2 * dad.k.D66 * (na1^2 - na2^2)
        g3a = (dad.k.D22 - dad.k.D12) * na1 * na2 + dad.k.D26 * (na1^2 - na2^2)


        g1d = (dad.k.D12 - dad.k.D11) * nd1 * nd2 + dad.k.D16 * (nd1^2 - nd2^2)
        g2d = 2 * (dad.k.D26 - dad.k.D16) * nd1 * nd2 + 2 * dad.k.D66 * (nd1^2 - nd2^2)
        g3d = (dad.k.D22 - dad.k.D12) * nd1 * nd2 + dad.k.D26 * (nd1^2 - nd2^2)


        d2tnadx2 = -(g1a * d4wdx4 + g2a * d4wdx3dy + g3a * d4wdx2dy2)
        d2tnddx2 = -(g1d * d4wdx4 + g2d * d4wdx3dy + g3d * d4wdx2dy2)

        d2tnady2 = -(g1a * d4wdx2dy2 + g2a * d4wdxdy3 + g3a * d4wdy4)
        d2tnddy2 = -(g1d * d4wdx2dy2 + g2d * d4wdxdy3 + g3d * d4wdy4)

        d2tnadxdy = -(g1a * d4wdx3dy + g2a * d4wdx2dy2 + g3a * d4wdxdy3)
        d2tnddxdy = -(g1d * d4wdx3dy + g2d * d4wdx2dy2 + g3d * d4wdxdy3)

        d2Rdx2 = (d2tnddx2 - d2tnadx2)
        d2Rdy2 = (d2tnddy2 - d2tnady2)
        d2Rdxdy = (d2tnddxdy - d2tnadxdy)


        d2R = [d2Rdx2 d2Rdy2 d2Rdxdy]'

        d2w = [d2wdx2 d2wdy2 d2wdxdy]'

        RS[:, i] = d2R
        WS[:, i] = d2w
    end
    return RS, WS
end

function compute_dtRw(pf, nf, dad, pre)
    m1 = nf[1]
    m2 = nf[2]

    t1 = -m2
    t2 = m1
    RS = zeros(size(dad.k.cantos, 1))
    WS = zeros(size(dad.k.cantos, 1))
    for ii = 1:size(dad.k.cantos, 1)
        pc = dad.k.cantos[ii, 2:3]
        na1 = dad.k.cantos[ii, 4]
        na2 = dad.k.cantos[ii, 5]
        nd1 = dad.k.cantos[ii, 6]
        nd2 = dad.k.cantos[ii, 7]
        rs = pc - pf
        #dad.k.Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]

        G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
        H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
        C1 =
            ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[1])
        C2 =
            ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[2])
        C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
        # a = constFS
        a = 1
        theta = atan(r2, r1)
        r = norm(rs)
        dRdx,
        dRdy,
        dSdx,
        dSdy,
        d2Rdx2,
        d2Rdxdy,
        d2Rdy2,
        d2Sdx2,
        d2Sdxdy,
        d2Sdy2,
        d3Rdx2dy,
        d3Rdx3,
        d3Rdxdy2,
        d3Rdy3,
        d3Sdx2dy,
        d3Sdx3,
        d3Sdxdy2,
        d3Sdy3,
        extra = pre

        for i = 1:2
            constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
            constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            dRdx[i] =
                2 *
                r *
                (cos(theta) + dad.k.d[i] * sin(theta)) *
                (
                    log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) - 2
                ) -
                4 *
                r *
                dad.k.e[i] *
                sin(theta) *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            dRdy[i] =
                2 *
                r *
                (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                    dad.k.e[i]^2 * sin(theta)
                ) *
                (
                    log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) - 2
                ) -
                4 *
                r *
                dad.k.e[i] *
                (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
                atan((dad.k.e[i] * sin(theta)), (cos(theta) + dad.k.d[i] * sin(theta)))

            dSdx[i] =
                r *
                dad.k.e[i] *
                sin(theta) *
                (
                    log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) - 2
                ) +
                2 *
                r *
                (cos(theta) + dad.k.d[i] * sin(theta)) *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            dSdy[i] =
                r *
                dad.k.e[i] *
                (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
                (
                    log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) - 2
                ) +
                2 *
                r *
                (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                    dad.k.e[i]^2 * sin(theta)
                ) *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            d2Rdx2[i] =
                2 * log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                )

            ####

            d2Rdxdy[i] =
                2 *
                dad.k.d[i] *
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) -
                4 *
                dad.k.e[i] *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
            ###
            d2Rdy2[i] =
                2 *
                (dad.k.d[i]^2 - dad.k.e[i]^2) *
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) -
                8 *
                dad.k.d[i] *
                dad.k.e[i] *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))


            d2Sdx2[i] =
                2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            d2Sdxdy[i] =
                dad.k.e[i] * (log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                )) +
                2 *
                dad.k.d[i] *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))


            d2Sdy2[i] =
                2 *
                dad.k.d[i] *
                dad.k.e[i] *
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) +
                2 *
                (dad.k.d[i]^2 - dad.k.e[i]^2) *
                atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            d3Rdx3[i] = 4 * (cos(theta) + dad.k.d[i] * sin(theta)) / (r * constA)

            # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdx2dy[i] =
                4 * (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (r * constA)

            # Cálculo da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdxdy2[i] =
                4 * (
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
                ) / (r * constA)

            # Cálculo da derivada terceira de R1 e R2 em relação a y conforme Shi e Bezine
            d3Rdy3[i] =
                4 * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) / (r * constA)

            # Cálculo da derivada terceira de S1 e S2 em relação a x conforme Shi e Bezine
            d3Sdx3[i] = -2 * dad.k.e[i] * sin(theta) / (r * constA)

            # Cálculo da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
            d3Sdx2dy[i] = 2 * dad.k.e[i] * cos(theta) / (r * constA)

            # Cálculo da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
            d3Sdxdy2[i] =
                2 *
                dad.k.e[i] *
                (
                    2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                ) / (r * constA)

            # Cálculo da derivada terceira de S1 e S2 em relação a y conforme Shi e Bezine
            d3Sdy3[i] =
                2 *
                dad.k.e[i] *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                ) / (r * constA)

        end

        dwdx = 1 / (8 * pi) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
        dwdy = 1 / (8 * pi) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))

        #dad.k.Derivadas segundas da solução fundamental
        d2wdx2 =
            1 / (8 * pi) * (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
        d2wdxdy =
            1 / (8 * pi) *
            (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
        d2wdy2 =
            1 / (8 * pi) * (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))

        #dad.k.Derivadas terceiras da solução fundamental
        d3wdx3 =
            1 / (8 * pi) * (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
        d3wdx2dy =
            1 / (8 * pi) *
            (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
        d3wdxdy2 =
            1 / (8 * pi) *
            (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
        d3wdy3 =
            1 / (8 * pi) * (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))


        g1a = (dad.k.D12 - dad.k.D11) * na1 * na2 + dad.k.D16 * (na1^2 - na2^2)
        g2a = 2 * (dad.k.D26 - dad.k.D16) * na1 * na2 + 2 * dad.k.D66 * (na1^2 - na2^2)
        g3a = (dad.k.D22 - dad.k.D12) * na1 * na2 + dad.k.D26 * (na1^2 - na2^2)


        g1d = (dad.k.D12 - dad.k.D11) * nd1 * nd2 + dad.k.D16 * (nd1^2 - nd2^2)
        g2d = 2 * (dad.k.D26 - dad.k.D16) * nd1 * nd2 + 2 * dad.k.D66 * (nd1^2 - nd2^2)
        g3d = (dad.k.D22 - dad.k.D12) * nd1 * nd2 + dad.k.D26 * (nd1^2 - nd2^2)

        dtndxa = -(g1a * d3wdx3 + g2a * d3wdx2dy + g3a * d3wdxdy2)
        dtndya = -(g1a * d3wdx2dy + g2a * d3wdxdy2 + g3a * d3wdy3)

        dtndxd = -(g1d * d3wdx3 + g2d * d3wdx2dy + g3d * d3wdxdy2)
        dtndyd = -(g1d * d3wdx2dy + g2d * d3wdxdy2 + g3d * d3wdy3)

        dtndta = dtndxa * t1 + dtndya * t2
        dtndtd = dtndxd * t1 + dtndyd * t2

        RS[ii] = -(dtndtd - dtndta) / dad.k.D22

        WS[ii] = -(dwdx * t1 + dwdy * t2) / dad.k.D22
        # @infiltrate
    end
    return RS, WS
end
function calc_d2w_cont(dad, forcas, desloc, dwdt_cont)
    n_fis = size(dad.NOS, 1)
    d2wnt = zeros(n_fis, 3)
    dwdt = zeros(3)
    normal_fonte = calc_normais(dad)

    for j = 1:length(dad.ELEM)#Laço dos elementos
        elem_j = dad.ELEM[j]#Laço dos elementos
        x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
        inds = [2(elem_j.indices) .- 1 2(elem_j.indices)]'[:]
        # Calcula os deslocamentos [2x1] no nó interno
        # w_elem = desloc[2(elem_j.indices).-1]
        dwdn_elem = desloc[2(elem_j.indices)]
        Mn_elem = forcas[2(elem_j.indices)]
        dwdt_elem = dwdt_cont[elem_j.indices]
        D = [
            dad.k.D11 dad.k.D12 dad.k.D16
            dad.k.D12 dad.k.D22 dad.k.D26
            dad.k.D16 dad.k.D26 dad.k.D66
        ]
        for i = 1:length(elem_j.ξs)
            eet = elem_j.ξs[i]
            N_geo, dN = calc_fforma(eet, elem_j)
            dxdqsi = dN' * x # dx/dξ & dy/dξ
            dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
            d2wdt2 = dN' * dwdt_elem / dgamadqsi
            d2wdtdn = dN' * dwdn_elem / dgamadqsi
            m = normal_fonte[elem_j.indices[i], 1]
            n = normal_fonte[elem_j.indices[i], 2]

            T = [
                m^2 n^2 2*m*n
                n^2 m^2 -2*m*n
                -m*n m*n m^2-n^2
            ]

            Dnt = T \ D / T'

            matD = [
                -Dnt[1, 1] -Dnt[1, 2] -2*Dnt[1, 3]
                -Dnt[1, 2] -Dnt[2, 2] -2*Dnt[2, 3]
                -Dnt[1, 3] -Dnt[2, 3] -2*Dnt[3, 3]
            ]
            S = inv(matD)
            A = [[-1 0 0]' S[:, 2] S[:, 3]]
            b =
                [-S[1, 1] * Mn_elem[i] d2wdt2 - S[1, 2] * Mn_elem[i] d2wdtdn -
                                                                     S[1, 3] * Mn_elem[i]]'
            y = inv(A) * b
            d2wdn2 = y[1]
            d2wnt[elem_j.indices[i], :] = [d2wdn2, d2wdt2, d2wdtdn]
            # @infiltrate
        end
    end
    d2wnt
end
function compute_tens_cont(d2wnt, dad)
    MATERIAL = dad.k.Material

    n_laminas = length(MATERIAL[:, 1])
    n_elem = length(dad.ELEM)
    nnos = nc(dad)
    thickness = MATERIAL[:, 6]
    z = zeros(n_laminas)
    tens_cont = zeros(n_elem, 9, n_laminas)
    sigma = zeros(nnos, 3, n_laminas)
    d2wdt = zeros(n_elem, 3, 3)
    StressXY = zeros(n_elem, 3, 3)

    MomentsXY = zeros(n_elem, 3, 3)
    Strain = zeros(3)
    z[1] = thickness[1]
    normal_fonte = calc_normais(dad)

    Qbar = Array{Float64}(undef, 3, 3)
    for lamina = 1:n_laminas
        if (lamina > 1)
            z[lamina] = z[lamina-1] + thickness[lamina]
        end
        E1 = MATERIAL[lamina, 2]
        E2 = MATERIAL[lamina, 3]
        G12 = MATERIAL[lamina, 4]
        ni12 = MATERIAL[lamina, 5]
        theta = MATERIAL[lamina, 7]
        theta = theta * pi / 180
        ni21 = ni12 * E2 / E1

        Q11 = E1 / (1 - ni12 * ni21)
        Q22 = E2 / (1 - ni12 * ni21)
        Q66 = G12
        Q16 = 0
        Q26 = 0
        Q12 = ni21 * E1 / (1 - ni12 * ni21)

        Q = [
            Q11 Q12 Q16
            Q12 Q22 Q26
            Q16 Q26 Q66
        ]

        m = cos(theta)
        n = sin(theta)

        T = [
            m^2 n^2 2*m*n
            n^2 m^2 -2*m*n
            -m*n m*n m^2-n^2
        ]

        Qbar = inv(T) * Q * inv(T')
        compliance = inv(Q)

    end

    for el = 1:n_elem
        elem = dad.ELEM[el]
        inds = elem.indices
        nn = normal_fonte[inds, :]
        for i = 1:3
            m = nn[i, 1]
            n = nn[i, 2]
            T = [
                m^2 n^2 2*m*n
                n^2 m^2 -2*m*n
                -m*n m*n m^2-n^2
            ]
            d2wdn2 = d2wnt[inds[i], 1]
            d2wdt2 = d2wnt[inds[i], 2]
            d2wdndt = d2wnt[inds[i], 3]
            d2wdnt = [
                d2wdn2 d2wdndt/2 0
                d2wdndt/2 d2wdt2 0
                0 0 0
            ]
            d2wdxy = T * d2wdnt * inv(T)
            d2w = [d2wdxy[1, 1] d2wdxy[2, 2] d2wdxy[1, 2]]
            MomentsXY[el, i, 1] =
                -(dad.k.D11 * d2w[1] + dad.k.D12 * d2w[2] + 2 * dad.k.D16 * d2w[3])
            MomentsXY[el, i, 2] =
                -(dad.k.D12 * d2w[1] + dad.k.D22 * d2w[2] + 2 * dad.k.D26 * d2w[3])
            MomentsXY[el, i, 3] =
                -(dad.k.D16 * d2w[1] + dad.k.D26 * d2w[2] + 2 * dad.k.D66 * d2w[3])
            Strain[1] = -z[1] * d2w[1]
            Strain[2] = -z[1] * d2w[2]
            Strain[3] = -z[1] * 2 * d2w[3]
            d2wdt[el, i, 1] = d2w[1]
            d2wdt[el, i, 2] = d2w[2]
            d2wdt[el, i, 3] = d2w[3]
            Stress = Qbar * Strain
            StressXY[el, i, 1] = Strain[1]
            StressXY[el, i, 2] = Strain[2]
            StressXY[el, i, 3] = Strain[3]



            for lamina = 1:n_laminas
                epsilon_n = z[lamina] * d2wdn2
                epsilon_t = z[lamina] * d2wdt2
                gamma_nt = z[lamina] * d2wdndt
                epnt = [
                    epsilon_n gamma_nt/2 0
                    gamma_nt/2 epsilon_t 0
                    0 0 0
                ]
                epxy = T * epnt * inv(T)
                sigmaxy = Qbar[lamina] * epxy
                tens_cont[el, 3*i-2:3*i, lamina] =
                    [sigmaxy[1, 1] sigmaxy[2, 2] sigmaxy[1, 2]]

            end
        end
    end

    # M = [MomentsXY[1:n_elem/4, :, 2]
    #     MomentsXY[1+n_elem/4:2*n_elem/4, :, 1]
    #     MomentsXY[1+2*n_elem/4:3*n_elem/4, :, 2]
    #     # MomentsXY[1+3*n_elem/4:end, :, 1]]

    for i = 1:n_elem
        elem = dad.ELEM[i]
        inds = elem.indices
        for j = 1:n_laminas
            sigma[inds[1], :, j] = tens_cont[i, 1:3, j]
            sigma[inds[2], :, j] = tens_cont[i, 4:6, j]
            sigma[inds[3], :, j] = tens_cont[i, 7:9, j]
        end
    end
    sigma, MomentsXY, d2wdt, StressXY
end


function compute_tens_int(d2w, dad)
    Material = dad.k.Material

    #Identifica o numero de laminas [de cada lado da llinha neutra]
    No_Laminas = length(Material[:, 1])
    No_PtosInternos = ni(dad)
    StrainXY = zeros(No_PtosInternos, 3, No_Laminas * 2)
    StressXY = zeros(No_PtosInternos, 3, No_Laminas * 2)
    StressLT = zeros(No_PtosInternos, 3, No_Laminas * 2)
    # Cálculo dos momentos
    MomentsXY = zeros(No_PtosInternos, 3)

    for i = 1:No_PtosInternos
        MomentsXY[i, 1] =
            -(dad.k.D11 * d2w[i, 1] + dad.k.D12 * d2w[i, 2] + 2 * dad.k.D16 * d2w[i, 3])
        MomentsXY[i, 2] =
            -(dad.k.D12 * d2w[i, 1] + dad.k.D22 * d2w[i, 2] + 2 * dad.k.D26 * d2w[i, 3])
        MomentsXY[i, 3] =
            -(dad.k.D16 * d2w[i, 1] + dad.k.D26 * d2w[i, 2] + 2 * dad.k.D66 * d2w[i, 3])
    end

    #Calcula as matrizes Q e Qbarra para cada lamina
    h = zeros(No_Laminas * 2)
    z = zeros(No_Laminas * 2)
    Matriz_Q = zeros(3, 3)
    Qbarra = zeros(3, 3)
    Matriz_T = zeros(3, 3)
    Strain = zeros(3)

    for i = 1:No_Laminas*2

        #Calcula as propriedades mecanicas da lamina

        if i <= No_Laminas
            k = No_Laminas - (i - 1)
        else
            k = i - No_Laminas
        end

        EL = Material[k, 2]
        ET = Material[k, 3]
        GLT = Material[k, 4]
        niLT = Material[k, 5]
        niTL = niLT * ET / EL
        h[i] = Material[k, 6]
        Tetha = Material[k, 7] * pi / 180

        #Calculo da cota z para as laminas acima e abaixo da linha neutra.

        if i <= No_Laminas
            if i == 1
                z[1] = 0
                for t = 1:No_Laminas
                    z[1] = z[1] - Material[t, 6]
                end
            else
                z[i] = z[i-1] + h[i-1]
            end
        else
            if i == No_Laminas + 1
                z[i] = h[i]
            else
                z[i] = z[i-1] + h[i]
            end
        end

        #Monta a matriz Q para cada lamina

        #Matriz_Q=[EL[i]/(1-niLT[i]*niTL[i])         ET[i]*niLT/(1-niLT[i]*niTL[i])   0
        #             ET[i]*niLT/(1-niLT[i]*niTL[i])    ET[i]/(1-niLT[i]*niTL[i])        0
        #             0                                 0                                GLT[i]]

        Matriz_Q[1, 1] = EL / (1 - niLT * niTL)
        Matriz_Q[2, 2] = ET / (1 - niLT * niTL)
        Matriz_Q[1, 2] = ET * niLT / (1 - niLT * niTL)
        Matriz_Q[2, 1] = Matriz_Q[1, 2]
        Matriz_Q[3, 3] = GLT
        Matriz_Q[1, 3] = 0
        Matriz_Q[2, 3] = 0
        Matriz_Q[3, 1] = 0
        Matriz_Q[3, 2] = 0

        #Monta a matriz T para cada lamina

        #Matriz_T = [ cos(Tetha[i])^2               sin(Tetha[i])^2                2*sin(Tetha[i])*cos(Tetha[i])
        #                sin(Tetha[i])^2               cos(Tetha[i])^2               -2*sin(Tetha[i])*cos(Tetha[i])
        #               -sin(Tetha[i])*cos(Tetha[i])   sin(Tetha[i])*cos(Tetha[i])    cos(Tetha[i])^2-sin(Tetha[i])^2]

        Matriz_T[1, 1] = cos(Tetha)^2
        Matriz_T[1, 2] = sin(Tetha)^2
        Matriz_T[1, 3] = 2 * sin(Tetha) * cos(Tetha)
        Matriz_T[2, 1] = sin(Tetha)^2
        Matriz_T[2, 2] = cos(Tetha)^2
        Matriz_T[2, 3] = -2 * sin(Tetha) * cos(Tetha)
        Matriz_T[3, 1] = -sin(Tetha) * cos(Tetha)
        Matriz_T[3, 2] = sin(Tetha) * cos(Tetha)
        Matriz_T[3, 3] = cos(Tetha)^2 - sin(Tetha)^2

        #Calcula a mattriz Qbarra para cada lamina

        #Matriz_Qbarra[i]= Matriz_invT[i]*Matriz_Q[i]*Matriz_T[i]

        Q11 = Matriz_Q[1, 1]
        Q12 = Matriz_Q[1, 2]
        Q22 = Matriz_Q[2, 2]
        Q66 = Matriz_Q[3, 3]

        Qbarra[1, 1] =
            Q11 * cos(Tetha)^4 +
            Q22 * sin(Tetha)^4 +
            2 * (Q12 + 2 * Q66) * sin(Tetha)^2 * cos(Tetha)^2
        Qbarra[2, 2] =
            Q11 * sin(Tetha)^4 +
            Q22 * cos(Tetha)^4 +
            2 * (Q12 + 2 * Q66) * sin(Tetha)^2 * cos(Tetha)^2
        Qbarra[1, 2] =
            (Q11 + Q22 - 4 * Q66) * sin(Tetha)^2 * cos(Tetha)^2 +
            Q12 * (cos(Tetha)^4 + sin(Tetha)^4)
        Qbarra[3, 3] =
            (Q11 + Q22 - 2 * Q12 - 2 * Q66) * sin(Tetha)^2 * cos(Tetha)^2 +
            Q66 * (sin(Tetha)^4 + cos(Tetha)^4)
        Qbarra[1, 3] =
            (Q11 - Q12 - 2 * Q66) * cos(Tetha)^3 * sin(Tetha) -
            (Q22 - Q12 - 2 * Q66) * cos(Tetha) * sin(Tetha)^3
        Qbarra[2, 3] =
            (Q11 - Q12 - 2 * Q66) * cos(Tetha) * sin(Tetha)^3 -
            (Q22 - Q12 - 2 * Q66) * cos(Tetha)^3 * sin(Tetha)
        Qbarra[2, 1] = Qbarra[1, 2]
        Qbarra[3, 2] = Qbarra[2, 3]
        Qbarra[3, 1] = Qbarra[1, 3]

        #Calcula as deformaçoes nos pontos internos

        for j = 1:No_PtosInternos
            Strain[1] = -z[i] * d2w[j, 1]
            Strain[2] = -z[i] * d2w[j, 2]
            Strain[3] = -z[i] * 2 * d2w[j, 3]

            #Calcula as tensoes nos pontos internos

            #Calcula as tensoes SigmaX; SigmaY e ThalXY nos pontos internos
            Matriz_StressXY = Qbarra * Strain

            #Calcula as tensoes SigmaL; SigmaT e ThalLT nos pontos internos
            Matriz_StressLT = Matriz_T * Matriz_StressXY


            # Cria o os vetores de saida

            StrainXY[j, 1:3, i] = Strain'
            StressXY[j, 1:3, i] = Matriz_StressXY'
            StressLT[j, 1:3, i] = Matriz_StressLT'
        end

    end

    MomentsXY, StrainXY, StressXY, StressLT
end

function compute_tens(d2w, d2wnt, dad)
    tens_cont, MomentsCont, d2w_cont, StressXY2 = compute_tens_cont(d2wnt, dad)
    MomentsXY, StrainXY, StressXY, StressLT = compute_tens_int(d2w, dad)
    # @infiltrate
    tens_cont, StressXY
end
