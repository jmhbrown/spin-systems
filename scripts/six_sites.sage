import os

cwd = os.getcwd()
print(cwd)
load('/home/jenny/code/spin-systems/interaction.sage')
load('/home/jenny/code/spin-systems/density_matrix.sage')

proj_a = projection(support_basis(density_matrix(3, t='A')))
proj_b = projection(support_basis(density_matrix(3, t='B')))

ham_A = hamiltonian(proj_a, length=6)
ham_B = hamiltonian(proj_b, length=6)

file_A = os.path.join(cwd, 'ham_A_s6.obj')
file_B = os.path.join(cwd, 'ham_B_s6.obj')

save(ham_A, file_A)
save(ham_B, file_B)

exit()

