cd "./build"
noise -g -f 0.7
L0min -v -l 0.01 -o ../examples/vertex_denoised_g7_e-2.obj
L0min -e -l 0.01 -o ../examples/edge_denoised_g7_e-2.obj

noise -n -f 0.7
L0min -v -l 0.01 -o ../examples/vertex_denoised_n7_e-2.obj
L0min -e -l 0.01 -o ../examples/edge_denoised_n7_e-2.obj

::L0min -v -l 0.001 -o ../examples/vertex_denoised.obj
::L0min -e -l 0.001 -o ../examples/edge_denoised.obj
pause