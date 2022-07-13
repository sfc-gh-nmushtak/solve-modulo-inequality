all_m_scales = []

# Appendix B of https://arxiv.org/pdf/2101.11408.pdf
# Generates all possible values of m_scale
for q in range(-342, -27):
    power5 = 5**-q
    z = 0
    while (1 << z) < power5:
        z += 1
    b = 2 * z + 2 * 64
    c = 2**b // power5 + 1
    while c >= (1 << 128):
        c //= 2
    all_m_scales.append(c)
for q in range(-27, 0):
    power5 = 5**-q
    z = 0
    while (1 << z) < power5:
        z += 1
    b = z + 127
    c = 2**b // power5 + 1
    all_m_scales.append(c)
for q in range(0, 308 + 1):
    power5 = 5**q
    while power5 < (1 << 127):
        power5 *= 2
    while power5 >= (1 << 128):
        power5 //= 2
    all_m_scales.append(power5)

# Returns the continued fraction of numer/denom as a list [a0; a1, a2, ..., an]
def continued_fraction(numer, denom):
    # If numer/denom is an integer, just return one term:
    if numer % denom == 0:
        return [numer // denom]
    # numer/denom = z + 1/f where z=floor(numer/denom)
    # -> f = 1/(numer/denom-z) = 1/((numer % denom)/denom) = denom/(numer % denom)
    return [numer // denom] + continued_fraction(denom, numer % denom)


# Given a continued fraction [a0; a1, a2, ..., an], returns the fraction equal to the continued fraction
def from_continued_fraction_to_fraction(cf):
    first_term = cf[0]
    # If the continued fraction has one term, then just return that term over 1:
    if len(cf) == 1:
        return first_term, 1

    r_numer, r_denom = from_continued_fraction_to_fraction(cf[1:])
    # numer/denom = first_term + 1/(r_numer/r_denom)
    #             = first_term + r_denom/r_numer
    #             = (r_numer*first_term + r_denom)/r_numer
    numer = r_numer * first_term + r_denom
    denom = r_numer

    # We don't need to reduce numer/denom: We can show gcd(numer, denom) = 1 by induction
    # If len(cf) == 1, then gcd(first_term, 1) == 1
    # Otherwise, gcd(numer, denom) = gcd(r_numer * first_term + r_denom, r_numer) = gcd(r_numer, r_denom) = 1 by induction hypothesis

    return numer, denom


# Given a continued fraction [a0; a1, a2, ..., an], returns all the convergents of that continued fraction
# as pairs of the form (numer, denom), where numer/denom is a convergent of the continued fraction
def convergents(cf):
    # If the continued fraction has one term, then just return that term over 1:
    if len(cf) == 1:
        return [(cf[0], 1)]
    # Otherwise, compute the convergents of [a0; a1, a2, ..., a(n-1)]
    # and then add the value of [a0; a1, a2, ..., an]
    return convergents(cf[:-1]) + [from_continued_fraction_to_fraction(cf)]


if __name__ == "__main__":
    # Enumerate through all the denominators of the convergents of m_scale / 2^137 up to 2^64
    for j, m_scale in enumerate(all_m_scales):
        for _, denom in convergents(continued_fraction(m_scale, 2**137)):
            if denom >= 2**64:
                break

            # Check if decimalSignificand=denom satisfies the inequality (decimalSignificand * m_scale) % 2^137 >= 2^137-2^64
            if (denom * m_scale) % (2**137) >= ((2**137) - (2**64)):
                print(
                    f"SOLUTION: decimalScale={j-342} m_scale={m_scale} decimalSignificand={denom}"
                )
