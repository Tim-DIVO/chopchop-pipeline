"""Global pairwise protein alignment, BLOSUM62 + affine gaps (Gotoh), stdlib only."""
B62="""A R N D C Q E G H I L K M F P S T W Y V B Z X *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1"""
_l=B62.strip().splitlines()
_h=_l[0].split()
M={}
for row in _l[1:]:
    p=row.split(); a=p[0]
    for b,v in zip(_h,p[1:]): M[(a,b)]=int(v)
def sc(a,b): return M.get((a,b), M.get((b,a),-1))

NEG = -10 ** 9


def align(s1: str, s2: str, go: int = -11, ge: int = -1):
    """Needleman-Wunsch with affine gaps (Gotoh). Returns the two gapped strings.

    go/ge are the BLAST/Clustal protein defaults (open -11, extend -1). For a
    pair of sequences this is the same computation Clustal Omega reduces to, so
    no external service and no e-mail address are needed."""
    n, m = len(s1), len(s2)
    Mx = [[NEG] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG] * (m + 1) for _ in range(n + 1)]   # gap in s2
    Iy = [[NEG] * (m + 1) for _ in range(n + 1)]   # gap in s1
    Mx[0][0] = 0
    for i in range(1, n + 1):
        Ix[i][0] = go + ge * (i - 1)
    for j in range(1, m + 1):
        Iy[0][j] = go + ge * (j - 1)
    for i in range(1, n + 1):
        a = s1[i - 1]
        Mi, Mp = Mx[i], Mx[i - 1]
        Xi, Xp = Ix[i], Ix[i - 1]
        Yi, Yp = Iy[i], Iy[i - 1]
        for j in range(1, m + 1):
            d = sc(a, s2[j - 1])
            Mi[j] = max(Mp[j - 1], Xp[j - 1], Yp[j - 1]) + d
            Xi[j] = max(Mp[j] + go, Xp[j] + ge, Yp[j] + go)
            Yi[j] = max(Mi[j - 1] + go, Yi[j - 1] + ge, Xi[j - 1] + go)
    i, j = n, m
    st = max((Mx[i][j], 'M'), (Ix[i][j], 'X'), (Iy[i][j], 'Y'))[1]
    A, B = [], []
    while i > 0 or j > 0:
        if i == 0:
            st = 'Y'
        elif j == 0:
            st = 'X'
        if st == 'M':
            A.append(s1[i - 1]); B.append(s2[j - 1])
            cur = Mx[i][j] - sc(s1[i - 1], s2[j - 1])
            i -= 1; j -= 1
            st = 'M' if cur == Mx[i][j] else ('X' if cur == Ix[i][j] else 'Y')
        elif st == 'X':
            A.append(s1[i - 1]); B.append('-')
            cur = Ix[i][j]; i -= 1
            st = 'M' if cur == Mx[i][j] + go else ('X' if cur == Ix[i][j] + ge else 'Y')
        else:
            A.append('-'); B.append(s2[j - 1])
            cur = Iy[i][j]; j -= 1
            st = 'M' if cur == Mx[i][j] + go else ('Y' if cur == Iy[i][j] + ge else 'X')
    return ''.join(reversed(A)), ''.join(reversed(B))


def residue_maps(a1: str, a2: str):
    """1-based residue index maps between the two ungapped sequences."""
    i = j = 0
    h2z: dict[int, int] = {}
    z2h: dict[int, int] = {}
    for c1, c2 in zip(a1, a2):
        if c1 != '-':
            i += 1
        if c2 != '-':
            j += 1
        if c1 != '-' and c2 != '-':
            h2z[i] = j
            z2h[j] = i
    return h2z, z2h


def conservation(a1: str, a2: str, zlen: int, win: int = 11):
    """Per-residue conservation of sequence 2 (the zebrafish one).

    For every residue of s2 we take a window of `win` alignment columns centred
    on its column and report (a) identity %, (b) BLOSUM-positive %, over the
    columns where both sequences have a residue. A residue with no aligned
    partner anywhere in its window gets None."""
    col_of = [None] * (zlen + 1)
    j = 0
    pair = []
    for n, (c1, c2) in enumerate(zip(a1, a2)):
        if c2 != '-':
            j += 1
            if j <= zlen:
                col_of[j] = n
        pair.append((c1, c2))
    half = win // 2
    ident, simil = [None] * (zlen + 1), [None] * (zlen + 1)
    for r in range(1, zlen + 1):
        c = col_of[r]
        if c is None:
            continue
        lo, hi = max(0, c - half), min(len(pair), c + half + 1)
        tot = idn = sim = 0
        for c1, c2 in pair[lo:hi]:
            if c1 == '-' or c2 == '-':
                continue
            tot += 1
            if c1 == c2:
                idn += 1
            if sc(c1, c2) > 0:
                sim += 1
        if tot:
            ident[r] = round(100.0 * idn / tot)
            simil[r] = round(100.0 * sim / tot)
    return ident, simil
