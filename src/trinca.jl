function encontra_trinca(MALHA, indica)
    n = size(MALHA, 1)
    conta = 1
    indices = []
    for i in 1:n
        contaf = conta + MALHA[i, 2]
        if indica[i] == true
            indices = [indices; conta:contaf-1]
        end
        conta = contaf
    end
    indices
end