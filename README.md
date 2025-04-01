# ML-Based Acoustic Feedback Path Estimation

## Overview
This repository presents an implementation of a **machine learning (ML)-based approach** for **acoustic feedback path estimation** in hearing aids. The project focuses on improving acoustic path modeling using **adaptive filtering techniques**, leveraging a **tap-amplitude-based partitioned block proportionate algorithm** alongside **decorrelation techniques** and **adaptive noise variance estimation**.

## Features
- **Partitioned Block Proportionate Adaptive Algorithm**:
  - Enhances acoustic path modeling for hearing aid applications.
  - Achieves **PESQ scores of 4.35 and 4.30**, demonstrating improved signal quality.
- **Decorrelation and Adaptive Noise Estimation**:
  - Implements a **real-time decorrelation algorithm**.
  - Integrates **online adaptive noise variance estimation**, resulting in a **20% improvement** in feedback path estimation.
- **Performance Analysis**:
  - Evaluates **transient and steady-state convergence** in the time domain.
  - Implements the algorithm in the **frequency domain** for comprehensive analysis.

## Implementation
- **Language**: MATLAB
- **Core Techniques**:
  - Frequency-Domain Adaptive Filtering
  - Partitioned Block Proportionate Algorithms
  - Real-Time Noise Variance Adaptation
  - Acoustic Path Model Optimization

## Usage
To run the project:
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/ml-acoustic-feedback.git
   cd ml-acoustic-feedback
   ```
2. Open MATLAB and navigate to the project directory.

*Codes used for evaluation purposes were not uploaded.
## Results
The proposed algorithm significantly enhances **acoustic feedback path estimation**, leading to superior performance in hearing aid applications. The results confirm **higher PESQ scores** and **improved noise variance modeling**, reducing feedback distortion.

## Future Work
- Implement real-time processing for embedded hearing aid systems.
- Optimize computational efficiency for low-power applications.
- Explore deep learning models for further refinement.

## Contributing
Contributions are welcome! Fork the repository and submit pull requests.

## License
Unethical use of the documents/codes can lead to legal action. Please contact the owner for collaboration or any other use. 

## Contact
For queries or collaborations, please reach out via [your email/contact information].

