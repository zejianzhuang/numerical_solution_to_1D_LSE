target:
	clear
	julia --project=. -e "@time include(\"./Dmat.jl\")"
