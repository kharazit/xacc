/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/
#ifndef OPENPULSEVISITOR_HPP_
#define OPENPULSEVISITOR_HPP_

#include <memory>
#include "AllGateVisitor.hpp"
#include "pulse_instruction.hpp"
#include "pulse_composite.hpp"
#include "json.hpp"

using nlohmann::json;

namespace xacc {
namespace quantum {

class OpenPulseVisitor: public BaseInstructionVisitor,
                       public InstructionVisitor<xacc::quantum::PulseInstruction>,
                       public InstructionVisitor<xacc::quantum::PulseComposite> {
protected:

	constexpr static double pi = 3.1415926;
    std::string instructions = "";

    int runningTime = 0;
    int currentTime = 0;

public:
    std::vector<json> instructionsJson;
    std::vector<json> pulseLibraryJson;

	const std::string name() const override {
		return "openpulse-visitor";
	}

	const std::string description() const override{
		return "Map XACC IR to OpenPulse.";
	}

	const std::string toString() override {
		return native;
	}


	void visit(PulseInstruction& i) override {
        std::vector<std::complex<double>> samples = i.getParameter(3).as<std::vector<std::complex<double>>>();
        std::vector<std::vector<double>> tmp;
        for (auto& s : samples) tmp.push_back({s.real(), s.imag()});
        json j, ij;
        j["name"] = i.name();
        if (i.name() == "fc") {
            j["phase"] = i.getParameter(2).as<double>();
        }
        // j["samples"] = tmp;
        pulseLibraryJson.push_back(j);

        ij["name"] = i.name();
        ij["ch"] = i.getParameter(0).toString();
        ij["phase"] = i.getParameter(2).as<double>();
        ij["t0"] = runningTime + i.getParameter(1).as<int>();
        instructionsJson.push_back(ij);

        currentTime = runningTime + i.getParameter(1).as<int>();
	}

    void visit(PulseComposite& c) override {
        runningTime += currentTime;
    }

	virtual ~OpenPulseVisitor() {}
};


}
}

#endif
