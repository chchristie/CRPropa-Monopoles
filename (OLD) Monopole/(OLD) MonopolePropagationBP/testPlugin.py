import crpropa
import MonopolePropagationBP

print("My Simulation\n")

ml = crpropa.ModuleList()

ml.add(MonopolePropagationBP.MonopolePropagationBP())
ml.add(crpropa.MaximumTrajectoryLength(1000*crpropa.parsec))

print("+++ List of modules")
print(ml.getDescription())


print("+++ Preparing source")
source = crpropa.Source()
print(source.getDescription())

print("+++ Starting Simulation")
ml.run(source, 1)

print("+++ Done")
