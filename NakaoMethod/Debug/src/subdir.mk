################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NakaoExperiment.cpp \
../src/NakaoExperiment2D.cpp \
../src/NakaoExperimentApprox.cpp \
../src/NakaoMethod.cpp \
../src/nintlib.cpp 

OBJS += \
./src/NakaoExperiment.o \
./src/NakaoExperiment2D.o \
./src/NakaoExperimentApprox.o \
./src/NakaoMethod.o \
./src/nintlib.o 

CPP_DEPS += \
./src/NakaoExperiment.d \
./src/NakaoExperiment2D.d \
./src/NakaoExperimentApprox.d \
./src/NakaoMethod.d \
./src/nintlib.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


