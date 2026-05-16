set(SOURCES
Chemistry.cxx
Configurator.cxx
Crc32.cxx
DecayChannel.cxx
DecayTable.cxx
Energy.cxx
Entropy.cxx
Event.cxx
EventGenerator.cxx
Integrator.cxx
Model.cxx
Model_BlastWave.cxx
Model_HRG.cxx
Parser.cxx
Particle.cxx
ParticleCoor.cxx
ParticleDB.cxx
ParticleDecayer.cxx
ParticleType.cxx
Pressure.cxx
SoundVelocity.cxx
Temperature.cxx
Thermodynamics.cxx
Vector3D.cxx
Viscosity.cxx)

add_library(therminator SHARED ${SOURCES})
target_link_libraries(therminator ${ROOT_LIBRARIES})


