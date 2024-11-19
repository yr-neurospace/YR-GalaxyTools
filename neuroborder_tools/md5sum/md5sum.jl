using ArgParse
using YRUtils

s = ArgParseSettings(description="Print or check MD5 checksums.")
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

    "--md5_file"
    help = "input MD5 file"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = false

    "--input_dir"
    help = "input directory (i.e. all files that require MD5 generation or verification must be within this directory)"
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

    "--quiet"
    help = "if given, do not exit with non-zero code when MD5 checking failed"
    action = "store_true"

    "--ignore_missing_files"
    help = "if given, files within MD5 file but not existed will be skipped"
    action = "store_true"
end

args = parse_args(s)
for k in keys(args)
    if typeof(args[k]) <: AbstractArray
        args[k] = YRUtils.BaseUtils.flatten_array(args[k])
    end
end
args["froms2tos"] = collect(zip(args["files"], args["names"]))

YRUtils.BaseUtils.file_transfer(first.(args["froms2tos"]), last.(args["froms2tos"]); ops=["soft", "hard", "copy"])

if isempty(args["md5_file"])
    # print MD5 checksums
    YRUtils.BaseUtils.md5_gen(last.(args["froms2tos"]), joinpath(args["output_dir"][1], "MD5.txt"))
else
    # check MD5 checksums
    YRUtils.BaseUtils.md5_check(args["md5_file"][1], joinpath(args["output_dir"][1], "md5_check.txt");
        ignore_missing_files=args["ignore_missing_files"], add_root_dir=false, quiet=args["quiet"])
end
