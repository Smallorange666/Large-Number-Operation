def process_sample_data(p, samples):
    results = []
    for a, b in samples:
        add = (a + b) % p
        sub = (a - b) % p
        mul = (a * b) % p
        inv = pow(a, -1, p)  # 使用 pow 函数计算模逆
        exp = pow(a, b, p)
        results.append((add, sub, mul, inv, exp))
    return results


def read_input_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    p = int(lines[1].strip())
    samples = []
    for i in range(2, len(lines), 3):
        a = int(lines[i].strip())
        b = int(lines[i+1].strip())
        samples.append((a, b))
    return p, samples


def write_output_file(filename, results):
    with open(filename, 'w') as f:
        for add, sub, mul, inv, exp in results:
            f.write(f"{add}\n{sub}\n{mul}\n{inv}\n{exp}\n\n")


def main():
    input_filename = "input.txt"
    output_filename = "output.txt"
    p, samples = read_input_file(input_filename)
    results = process_sample_data(p, samples)
    write_output_file(output_filename, results)


if __name__ == "__main__":
    main()
