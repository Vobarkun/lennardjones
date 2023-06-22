using GLMakie, LennardJones, ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--darktheme", "-d"
        action = :store_true
    "--borderless", "-b"
        action = :store_true
end
args = parse_args(s)

if args["darktheme"]
    LennardJones.darktheme!()
end

fig, plt, node, screen = LennardJones.main(decorated = !args["borderless"])
wait(screen)