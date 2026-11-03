import numpy as np

from src.envs.hirzebruch import count_c4


N = 39
S = (0, 1, 3, 7, 12, 22)


def build_incidence_matrix():
    data = np.zeros((N, N), dtype=np.uint8)
    for line in range(N):
        for s in S:
            point = (line + s) % N
            data[line, point] = 1
    return data


def main():
    data = build_incidence_matrix()

    row_degrees = data.sum(axis=1)
    col_degrees = data.sum(axis=0)
    c4 = count_c4(data)

    print("Sidon set:", S)
    print("row degrees:", sorted(set(map(int, row_degrees))))
    print("column degrees:", sorted(set(map(int, col_degrees))))
    print("edges:", int(data.sum()))
    print("C4s:", c4)

    if not np.all(row_degrees == 6):
        raise RuntimeError("row degrees are not all 6")
    if not np.all(col_degrees == 6):
        raise RuntimeError("column degrees are not all 6")
    if c4 != 0:
        raise RuntimeError("cyclic construction is not C4-free")

    np.save("hirzebruch39_cyclic.npy", data)
    np.savetxt("hirzebruch39_cyclic.txt", data, fmt="%d")
    print("saved hirzebruch39_cyclic.npy and hirzebruch39_cyclic.txt")


if __name__ == "__main__":
    main()
