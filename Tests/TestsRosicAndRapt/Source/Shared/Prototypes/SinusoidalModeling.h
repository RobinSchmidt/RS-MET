#pragma once


template<class T>
std::vector<T> synthesizeSinusoidal(const RAPT::rsSinusoidalModel<T>& model, T sampleRate);


template<class T>
RAPT::rsSinusoidalModel<T> analyzeSinusoidal(T* sampleData, int numSamples, T sampleRate);