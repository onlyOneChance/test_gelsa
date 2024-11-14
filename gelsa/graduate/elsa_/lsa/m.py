import subprocess
import time

def main():
    script_name = "lsa_compute"
    input_file1 = "first_file.txt"
    input_file2 = "second_file.txt"
    result_file = "result"
    d_value = "10"
    r_value = "1"
    s_value = "50"
    fillMethod = 'quadratic'
    transFunc = 'simple'
    normMethod = 'percentile'
    pvalueMethod = ["theo","ddlsa","perm","bblsa","stlta"]
    # pvalueMethod = ["stlta"]
    precision = '1000'
    minOccur = '50'
    bootNum = '200'
    qvalueMethod = 'scipy'
    trendThresh = "0.1"
    approxVar = '1'

    for pm in pvalueMethod:
        print("###############################")
        print(f"{pm}_way计算p值")

        command = f"{script_name} {input_file1} {result_file} -e {input_file2} \
        -d {d_value} -r {r_value} -s {s_value} -p {pm} -T {trendThresh}"
        # -f {fillMethod} -t {transFunc} -n {normMethod} \
        # -m {minOccur} -x {precision} -b {bootNum} -q {qvalueMethod} \
        # -T {trendThresh} -a {approxVar} -v {progressive}"
        subprocess.run(command, shell=True)
        print(f"##########   {pm}    ##########")

    with open("result", "r") as f_in, open("result.lsa", "w") as f_out:
        for line in f_in:
            data = line.strip().split(",")
            formatted_line = ""
            for item in data:
                formatted_line += "{:<30}".format(item)
            f_out.write(formatted_line + "\n")


if __name__=="__main__":
    main()
# pip install numpy scipy pandas statsmodels 
