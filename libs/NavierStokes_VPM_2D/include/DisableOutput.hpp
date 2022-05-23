#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace VPM {

class DisableOutput {
public:
    DisableOutput(std::ostream& stream_) : stream(stream_) {
        oldbuf = stream.rdbuf();
        stream.rdbuf(temporaryOutput.rdbuf());
    }

    ~DisableOutput() {
        stream.rdbuf(oldbuf);
    }

private:
    std::ostream& stream;
    // Now we simply store it temporarily to a stringstream
    // TODO: Make a nullstream that ignores all output and does not store it.
    std::stringstream temporaryOutput;
    std::streambuf* oldbuf;
};
}
