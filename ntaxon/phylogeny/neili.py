import numpy as np
from scipy.spatial.distance import squareform

def _restdist_sitesort(sites, y, spp, alias, aliasweight):
    # Shell sort keeping alias, aliasweight in same order
    gap = int(sites / 2)
    while gap > 0:
        for i in range(gap + 1, sites + 1):
            j = int(i - gap)
            flip = True
            while (j > 0) & flip:
                jj = int(alias[j])
                jg = int(alias[j + gap])
                flip = False
                tied = True
                k = 1
                while (k <= spp) & tied:
                    flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1])
                    tied = (tied & (y[k - 1][jj - 1] == y[k - 1][jg - 1]))
                    k += 1

                if tied:
                    aliasweight[j] += aliasweight[j + gap]
                    aliasweight[j + gap] = 0

                if ~flip:
                    break

                itemp = alias[j]
                alias[j] = alias[j + gap]
                alias[j + gap] = itemp
                itemp = aliasweight[j]
                aliasweight[j] = aliasweight[j + gap]
                aliasweight[j + gap] = itemp
                j -= gap
        gap = int(gap / 2)

    return alias.astype("int"), aliasweight.astype("int")


def _restdist_sitecombine(sites, y, spp, alias, aliasweight):
    # combine sites that have identical patterns
    i = 1
    while i < sites:
        j = i + 1
        tied = True
        while (j <= sites) & tied:
            k = 1
            while (k <= spp) & tied:
                tied = (tied & (y[k - 1][alias[i] - 1] == y[k - 1][alias[j] - 1]))
                k += 1
            if tied & (aliasweight[j] > 0):
                aliasweight[i] += aliasweight[j]
                aliasweight[j] = 0
                alias[j] = alias[i]
            j += 1
        i = j - 1
    return alias.astype('int'), aliasweight.astype('int')

# TODO: Implement inputweights2 function properly for nei li distance
def _inputweights2(a, b, weight, prog):
    # /* input the character weights, 0 or 1 */
    ch = str()
    # long i;
    weightsum = 0

    for i in range(a, b):
        # do {
        #    if (eoln(weightfile))
        #    scan_eoln(weightfile);
        #    ch = gettc(weightfile);
        # } while (ch == ' ');

        weight[i] = 1

        if ch == '0' | ch == '1':
            weight[i] = ch - '0'
        else:
            print(f"ERROR: Bad weight character: {ch} -- weights in {prog} must be 0 or 1")
        weightsum += weight[i]
    # weights = True
    return True, weightsum  # weights, weightsum


def sitescrunch2(sites, i, j, alias, aliasweight):
    # move so positively weighted sites come first */

    done = False

    while done is False:
        if aliasweight[i - 1] > 0:
            i += 1
        else:
            if j <= i:
                j = i + 1
            if j <= sites:
                found = bool()
                # do {found = (aliasweight[j - 1] > 0);j++;} while (!(found || j > sites));
                while True:
                    found = (aliasweight[j - 1] > 0)
                    j += 1
                    if found | (j > sites):
                        break
                if found:
                    j -= 1
                    itemp = alias[i - 1]
                    alias[i - 1] = alias[j - 1]
                    alias[j - 1] = itemp
                    itemp = aliasweight[i - 1]
                    aliasweight[i - 1] = aliasweight[j - 1]
                    aliasweight[j - 1] = itemp
                else:
                    done = True
            else:
                done = True
        done = done | (i >= sites)

    return alias.astype('int'), aliasweight.astype('int')


def neili(X, neili=False, restsites=True, sitelength=6, ttratio=2.0, gama=False, cvi=0.0):
    """
    Calculates ``Nei and Li (1979)`` genetic distnace on a binary matrix (preferably from RFLP test),
    Containing sample as Rows and haplotypes as Columns

    Parameters
    ----------
    X : array_like
        An m by n array of m original observations in an
        n-dimensional space, containing value 1 for presence of an allele and 0 for absence

    restsites : bool, optional
        Restriction sites or fragments, default is True assuming binary data is for Restriction Sites.

    neili : bool, optional
        Original or modified Nei/Li model. True for Original Nei and Li (1979), and False for Modified
        Default is set for Modified Nei/Li

    sitelength : int, optional
        Restriction site length

    ttratio: float, optional
        Transition/transversion ratio - a real number greater than 0.0, as the expected ratio of transitions
        to transversions. Note that this is the resulting expected ratio of transitions to transversions.
        The default value is 2.0. Theis option is not available when you choose the original Nei/Li restriction
        fragment distance, which assumes a Jukes-Cantor (1969) model of DNA change, for which
        the transition/transversion ratio is in effect fixed at 0.5.

    gama: bool, optional
        Gamma distribution of rates among sites. This is different from the parameters used by Nei and Jin,
        who introduced Gamma distribution of rates in DNA distances, but related to their parameters:
        their parameter a is also known as "alpha", the shape parameter of the Gamma distribution.
        It is related to the coefficient of variation by ``CV = 1/(a^(1/2))`` or ``a = 1/(CV)^2``.
        (their parameter b is absorbed here by the requirement that time is scaled so that the mean rate
        of evolution is 1 per unit time, which means that a = b).
        As we consider cases in which the rates are less variable we should set a larger and larger,
        as CV gets smaller and smaller.
        The Gamma distribution option is not available when using the original Nei/Li restriction fragments distance.

    cvi: float, optional
        Coefficient of variation of substitution rate among sites (must be positive)

    Returns
    -------
    Y : ndarray
        Returns a condensed distance matrix Y.

    """

    y = np.ascontiguousarray(X)

    # if type(X) != np.ndarray:
    #    y = np.array(X)

    initialv = 0.1  # starting value  of branch length
    iterationsr = 20  # how many Newton iterations per distance

    spp, sites = y.shape  # species count, sites count

    if sitelength < 1.0:
        raise Exception("BAD RESTRICTION SITE LENGTH")



    # Input Options weights
    # TODO: Implement weights and weightsum
    weights = False
    weightsum = int()

    # Options for gamma
    if gama:
        if cvi <= 0:
            raise Exception("Coefficient of variation of substitution rate among sites")

        cvi = 1.0 / (cvi * cvi)

    xi = (ttratio - 0.5) / (ttratio + 0.5)
    xv = 1.0 - xi
    fracchange = xi * 0.5 + xv * 0.75

    # Assign and declare weights
    weight = np.zeros(sites + 1)
    alias = np.zeros(sites + 1)
    aliasweight = np.zeros(sites + 1)

    endsite = 0

    # Distance Matrix
    d = np.zeros((spp, spp))

    # From Inputoptions in phylip
    for i in range(1, sites + 1):
        weight[i] = 1

    weightsum = sites

    # make up weights vector to avoid duplicate computations
    for i in range(1, sites + 1):
        alias[i] = i
        aliasweight[i] = weight[i]

    # Reassign alias, aliasweight
    alias, aliasweight = _restdist_sitesort(sites, y, spp, alias, aliasweight)

    alias, aliasweight = _restdist_sitecombine(sites, y, spp, alias, aliasweight)
    alias, aliasweight = sitescrunch2(sites + 1, 2, 3, alias, aliasweight)

    # Makeweights
    for i in range(1, sites + 1):
        weight[i] = aliasweight[i]
        if weight[i] > 0:
            endsite = i
    weight[0] = 1

    def makev(m, n):
        # compute one distance

        g = 0
        # vv = 0
        numerator = 0
        denominator = 0
        for i in range(0, endsite):
            ii = alias[i + 1]
            if (y[m - 1][ii - 1] == 1) | (y[n - 1][ii - 1] == 1):
                denominator += weight[i + 1]
                if (y[m - 1][ii - 1] == 1) & (y[n - 1][ii - 1] == 1):
                    numerator += weight[i + 1]

        f = 2 * numerator / (denominator + numerator)
        if restsites:
            if np.exp(-1 * sitelength * 1.38629436) > f:
                print(f"ERROR: Infinite distance between species {m} and {n}")
                raise Exception(f"ERROR: Infinite distance between species {m} and {n}")

        if not restsites:
            if not neili:
                f = (np.sqrt(f * (f + 8.0)) - f) / 2.0
            else:
                g = initialv
                delta = g
                it = 0
                while np.fabs(delta) > 0.00002 and it < iterationsr:
                    it = it + 1  # it++
                    h = g
                    g = np.exp(0.25 * np.log(f * (3 - 2 * g)))
                    delta = g - h

        if restsites == False and neili:
            vv = - (2.0 / sitelength) * np.log(g)
        else:
            if neili == True and restsites == True:
                pp = np.exp(np.log(f) / (2 * sitelength))
                vv = -(3.0 / 2.0) * np.log((4.0 / 3.0) * pp - (1.0 / 3.0))
            else:
                pp = np.exp(np.log(f) / sitelength)
                delta = initialv
                tt = delta
                it = 0
                while np.fabs(delta) > 0.000001 and it < iterationsr:
                    it = it + 1  # it++;
                    if gama == True:
                        p1 = np.exp(-1 * cvi * np.log(1 + tt / cvi))
                        p2 = np.exp(-1 * cvi * np.log(1 + xv * tt / cvi)) - np.exp(-1 * cvi * np.log(1 + tt / cvi))
                        p3 = 1.0 - np.exp(-1 * cvi * np.log(1 + xv * tt / cvi))
                    else:
                        p1 = np.exp(-tt)
                        p2 = np.exp(-xv * tt) - np.exp(-tt)
                        p3 = 1.0 - np.exp(-xv * tt)
                    q1 = p1 + p2 / 2.0 + p3 / 4.0
                    g = q1 - pp
                    if g < 0:
                        delta = np.fabs(delta) / -2.0
                    else:
                        delta = np.fabs(delta)
                    tt += delta
                vv = fracchange * tt
        v = np.fabs(vv)
        return v

    # compute distance matrix
    for i in range(1, spp + 1):
        for j in range(i + 1, spp + 1):
            v = makev(i, j)
            d[i - 1][j - 1] = v
            d[j - 1][i - 1] = v

    return squareform(d)
