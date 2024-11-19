using ArgParse
using YRUtils

s = ArgParseSettings(description="Run Trim Galore on FASTQ files.")
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

dict = YRUtils.BioUtils.auto_detect_fastq_read_type(last.(args["froms2tos"]); allow_mix=false)
files_dict = if dict["paired"]["status"] == "yes"
    dict["paired"]["dict"]
elseif dict["single"]["status"] == "yes"
    dict["single"]["dict"]
else
    @error "did not detect any paired-end or single-end files"
end
files_read_type = if dict["paired"]["status"] == "yes"
    "paired"
elseif dict["single"]["status"] == "yes"
    "single"
else
    @error "did not detect any paired-end or single-end files"
end

YRUtils.BioUtils.trimgalore(files_dict, files_read_type, args["output_dir"][1];
    trimgalore_options="--cores 4 --phred33 --quality 20 --length 30 --trim-n",
    num_jobs=1)
