import random
from sympy import nextprime

def generate_prime(bits):
    # 生成一个随机的初始值
    start = random.getrandbits(bits)
    # 找到大于等于 start 的下一个素数
    prime = nextprime(start)
    return prime

def generate_sample_data(p, n):
    samples = []
    for _ in range(n):
        a = random.randint(1, p-1)
        b = random.randint(1, p-1)
        samples.append((a, b))
    return samples

def write_to_file(filename, p, samples):
    with open(filename, 'w') as f:
        f.write(f"{len(samples)}\n")
        f.write(f"{p}\n")
        for a, b in samples:
            f.write(f"{a}\n{b}\n\n")

def main(P_BIT, N):
    p = generate_prime(P_BIT)
    samples = generate_sample_data(p, N)
    write_to_file("input.txt", p, samples)

if __name__ == "__main__":
    P_BIT = 2048  # 可以根据需要调整
    N = 10  # 可以根据需要调整
    main(P_BIT, N)