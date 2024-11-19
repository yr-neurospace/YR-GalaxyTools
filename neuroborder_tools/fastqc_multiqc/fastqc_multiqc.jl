using ArgParse
using YRUtils

s = ArgParseSettings(description="Run fastqc/multiqc on FASTQ files.")
@add_arg_table s begin
    "--names"
    help = "input file names separated by white space"
    nargs = '+'
    arg_type = String
    action = "append_arg"
    required = true

    "--files"
    help = "input file paths separated by white space"
    nargs = '+'
    arg_type = String
    action = "append_arg"
    required = true

    "--input_dir"
    help = "input directory (i.e. all FASTQ files must be within this directory)"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--output_dir"
    help = "output directory"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true
end

args = parse_args(s)
for k in keys(args)
    if typeof(args[k]) <: AbstractArray
        args[k] = YRUtils.BaseUtils.flatten_array(args[k])
    end
end
args["froms2tos"] = collect(zip(args["files"], joinpath.(args["input_dir"][1], args["names"])))

YRUtils.BaseUtils.file_transfer(first.(args["froms2tos"]), last.(args["froms2tos"]); ops=["soft", "hard", "copy"])

YRUtils.BioUtils.fastqc(last.(args["froms2tos"]), args["output_dir"][1];
    fastqc_options="--threads 4", multiqc_options="--zip-data-dir", num_jobs=4)
