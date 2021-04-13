def inputweights2(a, b, weight, prog):
    #/* input the character weights, 0 or 1 */
    ch = str()
    #long i;
    weightsum = 0

    for i in range(a, b):
        #do {
        #    if (eoln(weightfile))
        #    scan_eoln(weightfile);
        #    ch = gettc(weightfile);
        #} while (ch == ' ');

        weight[i] = 1

        if ch == '0' | ch == '1':
            weight[i] = ch - '0'
        else:
            print(f"ERROR: Bad weight character: {ch} -- weights in {prog} must be 0 or 1")
        weightsum += weight[i]
    # weights = True
    return True, weightsum # weights, weightsum

def sitescrunch2(sites, i, j, alias, aliasweight):
    # move so positively weighted sites come first */

    done = False

    while ~done:
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
