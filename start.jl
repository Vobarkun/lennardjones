using GLMakie, LennardJones

GLMakie.activate!(renderloop = GLMakie.renderloop, framerate = 60, fxaa = false, ssao = false);
LennardJones.main();