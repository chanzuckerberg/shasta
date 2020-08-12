#ifdef SHASTA_HTTP_SERVER

// Implementation of class HttpServer - see HttpServer.hpp for more information.

// Shasta.
#include "HttpServer.hpp"
#include "filesystem.hpp"
#include "platformDependent.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ip/v6_only.hpp>
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
using namespace boost;
using namespace asio;
using namespace ip;

// Standard library.
#include "chrono.hpp"
#include "fstream.hpp"
#include "iostream.hpp"
#include <regex>
#include <sstream>
#include "stdexcept.hpp"

// Operating system.
#include <sys/types.h>
#include <unistd.h>


// This function puts the server into an endless loop
// of processing requests.
// This is the function that the base class should call to start the server.
void HttpServer::explore(uint16_t port, bool localOnly, bool sameUserOnly)
{
    // Sanity check on the arguments.
    if(!localOnly && sameUserOnly) {
        throw runtime_error("Http::explore called with localOnly=false and sameUserOnly=true. "
            "This combination is not allowed.");
    }

    // Create the acceptor, making sure to accept both ipv4 and ipv6 ip addresses.
    io_service service;
    tcp::acceptor acceptor(service);
    tcp::endpoint endpoint = (
        localOnly ?
        tcp::endpoint(ip::address::from_string("::ffff:127.0.0.1"), port) :
        tcp::endpoint(tcp::v6(), port)
    );
    acceptor.open(endpoint.protocol());
    v6_only ipv6Option(false);
    acceptor.set_option(ipv6Option);

    // Bind to the requested port, and try the next port if that fails.
    bool bindWasSuccessful = false;
    for(int iteration=0; iteration<30; ++iteration) {
        try {
            acceptor.bind(endpoint);
            bindWasSuccessful = true;
            break;
        } catch(...) {
            // The bind failed. try again on the next port.
            cout << "Port " << port << " is not available." << endl;
            endpoint.port(++port);
        }
    }
    if(!bindWasSuccessful) {
        throw runtime_error("Unable to find a usable port.");
    }

    // The acceptor is bound to this port. Start listening for connections.
    acceptor.listen();
    cout << "Listening for http requests on port " << port << endl;
    if(localOnly) {
        cout << "To connect, "
            "point your browser to http://localhost:" << port << endl;
        if(sameUserOnly) {
            cout << "Only accepting local connections originating from a process "
                "owned by the same user running the server." << endl;
            // Also start the default browser and point it to the server.
#ifdef __linux__
            ::system(("xdg-open http://localhost:" + to_string(port)).c_str());
#else
            ::system(("open http://localhost:" + to_string(port)).c_str());
#endif
        } else {
            cout << "Accepting local connections from any user. "
                "Connections from the local computer are accepted from any user. "
                "This means that other users with access to the local computer can "
                "access the server. "
                "THIS CHOICE ALLOWS OTHER USERS TO LOOK AT YOUR DATA AND SHOULD "
                "NOT BE USED IF YOUR ASSEMBLY DATA IS SUBJECT TO "
                "CONFIDENTIALITY RESTRICTIONS OR IS NOT CLEARED OR CONSENTED "
                "FOR PUBLIC RELEASE."
                << endl;
        }
    } else {
        cout << "Accepting connections from any computer, any user. " << endl;
        cout << "To connect from this computer, "
            "point your browser to http://localhost:" << port << endl;
        cout << "To connect from a different computer, "
            "point your browser to http://XYZ:" << port << 
            ", where XYZ is either the name or IP address of this computer." << endl;
        cout << 
            "All connections are accepted: local and remote, from any user. "
            "This means that all users, not only on the computer running the server, "
            "but also on all computers on the same local area network, can use the server. "
            "In addition, if the computer running the server is not protected by a firewall, "
            "everybody on the Internet can also access the server. "
            "THIS CHOICE ALLOWS OTHER USERS, EVEN ON DIFFERENT COMPUTERS, "
            "TO LOOK AT YOUR DATA AND SHOULD "
            "NOT BE USED IF YOUR ASSEMBLY DATA IS SUBJECT TO "
            "CONFIDENTIALITY RESTRICTIONS OR IS NOT CLEARED OR CONSENTED "
            "FOR PUBLIC RELEASE." 
            << endl;
    }



    // Endless loop over incoming connections.
    while(true) {
          tcp::iostream s;
          tcp::endpoint remoteEndpoint;
          boost::system::error_code errorCode;
          acceptor.accept(*s.rdbuf(), remoteEndpoint, errorCode);
          if(errorCode) {
              // If interrupted with Ctrl-C, we get here.
              cout << "\nError code from accept: " << errorCode.message() << endl;
              s.close();        // Should not be necessary.
              acceptor.close(); // Should not be necessary
              return;
          }

          // If sameUserOnly was specified, check that this is a local
          // connection originating from a process owned by the same
          // user running the server.
          if(sameUserOnly) {
              SHASTA_ASSERT(localOnly);
              if(!isLocalConnectionSameUser(s, port)) {
                  // Unceremoniously close the connection.
                  cout << timestamp << "Reset a local connection originating from a process "
                      "not owned by the same user running the server." << endl;
                  continue;
              }
          }

          // Process the request.
          cout << timestamp << remoteEndpoint.address().to_string() << " " << flush;
          const auto t0 = steady_clock::now();
          processRequest(s);
          const auto t1 = steady_clock::now();
          cout << timestamp << "Request satisfied in " << seconds(t1 - t0) << "s." << endl;
    }
}


void HttpServer::setRequestTimeout(int tsec, tcp::iostream& s) {
#if BOOST_VERSION < 106600
    s.expires_from_now(boost::posix_time::seconds(tsec));
#else
    s.expires_after(std::chrono::seconds(tsec));
#endif
}

void HttpServer::processRequest(tcp::iostream& s)
{
    // If the client is too slow sending the request, drop it.
    setRequestTimeout(1, s);

    // Get the first line, which must contain the GET request.
    string requestLine;
    getline(s, requestLine);
    if(requestLine.empty()) {
        cout << "Empty request ignored." << endl;
        return;
    }

    // Parse it to get only the request string portion.
    // It is the second word of the first line.
    vector<string> tokens;
    boost::algorithm::split(tokens, requestLine, boost::algorithm::is_any_of(" "));
    if(tokens.size() != 3) {
        s << "Unexpected number of tokens in http request: expected 3, got " << tokens.size();
        cout << "Unexpected number of tokens in http request: expected 3, got " << tokens.size() << endl;
        cout << "Request was: " << requestLine << endl;
        return;
    }
    if(tokens.front() == "POST") {
        setRequestTimeout(10000000, s);
        processPost(tokens, s);
        return;
    }

    if(tokens.front() != "GET") {
        s << "Unexpected keyword in http request: " << tokens.front();
        cout << "Unexpected keyword in http request: " << tokens.front() << endl;
        cout << "Request was: " << requestLine << endl;
        return;
    }
    const string& request = tokens[1];
    if(request.empty()) {
        s << "Empty GET request: " << requestLine;
        cout << "Empty GET request: " << requestLine;
        return;
    }

    // Give ourselves time to satisfy the request
    setRequestTimeout(86400, s);

    // Parse the request.
    cout << requestLine << endl;
    boost::algorithm::split(tokens, request, boost::algorithm::is_any_of("?=&"));

    // Do URl decoding on each token.
    // This takes care of % encoding, which the browser will do if it has to send special characters.
    // Note that we have to do this after parsing the request into tokens.
    // With this, we can support special characters in cell meta data, cell set names, graph names, etc.
    for(string& token : tokens) {
        string newToken;
        urlDecode(token, newToken);
        if(newToken != token) {
            cout << "Request token " << token << " decoded as " << newToken << endl;
        }
        token = newToken;
    }



    // Read the rest of the input from the client, but ignore it,
    // except for the User Agent string, which tells us what browser
    // issued the request.
    // If we don't read all the input, the client may get a timeout.
    string line;
    const string userAgentPrefix = "User-Agent: ";
    BrowserInformation browserInformation;
    string originatingProcessUserName;
    string thisProcessUserName;
    while(true) {
        if(!s) {
            break;
        }
        getline(s, line);
        if(!s) {
            break;
        }
        // cout << "Got another line of length " << line.size() << ": <<<" << line << ">>>" << endl;
        if(line.size()==1) {
            break;
        }

        // See if this is the User Agent string.
        if(line.compare(0, userAgentPrefix.size(), userAgentPrefix) == 0) {
            browserInformation.set(line);
        }
    }
    cout << "isFirefox=" << browserInformation.isFirefox << " ";
    cout << "isChrome=" << browserInformation.isChrome << endl;



    // Write the success response.
    // We don't write the required empty line, so the derived class can send headers
    // if it wants to.
    s << "HTTP/1.1 200 OK\r\n";

    // The derived class processes the request.
    processRequest(tokens, s, browserInformation);
}



#if 0
// If the derived class wants to process POST requests,
// it should override this. The version implemented in this class
// just returns an error.
void HttpServer::processPost(
    const vector<string>& requestLine,
    const string& postData,
    std::ostream& s )
{
    // s << "HTTP/1.1 405 Method Not Allowed\r\n\r\n";
    cout << "***abc" << endl;
    s << "HTTP/1.1 200 OK\r\n\r\n";

    // Copy to a file.
    ofstream debugOut("Post.txt");
    debugOut << postData;



}
#endif



// Construct browser information from the User Agent line.
void HttpServer::BrowserInformation::set(const string& userAgentHeader)
{
    const string firefoxString = "Firefox/";
    const string chromeString = "Chrome/";
    const string edgeString = "Edge/";

    // Tokenize the user agent header.
    boost::tokenizer< boost::char_separator<char> > tokenizer(userAgentHeader, boost::char_separator<char>(" "));
    vector<string> tokens;
    tokens.insert(tokens.begin(), tokenizer.begin(), tokenizer.end());

    for(const string& token: tokens) {
        if(token.compare(0, firefoxString.size(), firefoxString) == 0) {
            isFirefox = true;
        }
        if(token.compare(0, chromeString.size(), chromeString) == 0) {
            isChrome = true;
        }
        if(token.compare(0, edgeString.size(), edgeString) == 0) {
            isEdge = true;
        }
    }

    // Edge writes both Chrome and Edge in its User Agent string.
    if(isChrome && isEdge) {
        isChrome = false;
    }

}



// Return all values assigned to a parameter.
// For example, if the request has ...&a=xyz&a=uv,
// when called with argument "a" returns a set containing "xyz" and "uv".
void HttpServer::getParameterValues(
    const vector<string>& request,
    const string& name,
    std::set<string>& values)
{
    for(size_t i=0; i<request.size()-1; i++) {
        if(request[i]==name) {
            values.insert(request[i+1]);
        }
    }

}
void HttpServer::getParameterValues(const vector<string>& request, const string& name, vector<string>& values)
{
    for(size_t i=0; i<request.size()-1; i++) {
        if(request[i]==name) {
            values.push_back(request[i+1]);
        }
    }

}


// This takes care of percent encoding in the url.
// I adapted it from boost asio example http/server3.
bool HttpServer::urlDecode(const string& in, string& out)
{
  out.clear();
  out.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i)
  {
    if (in[i] == '%')
    {
      if (i + 3 <= in.size())
      {
        int value = 0;
        std::istringstream is(in.substr(i + 1, 2));
        if (is >> hex >> value)
        {
          out += static_cast<char>(value);
          i += 2;
        }
        else
        {
          return false;
        }
      }
      else
      {
        return false;
      }
    }
    else if (in[i] == '+')
    {
      out += ' ';
    }
    else
    {
      out += in[i];
    }
  }
  return true;
}



// Function to do URL encoding.
// This is necessary when we create a hyperlink that contains special characters.
// Adapted from here:
// https://stackoverflow.com/questions/154536/encode-decode-urls-in-c
string HttpServer::urlEncode(const string& s)
{
    std::ostringstream escaped;
    escaped.fill('0');
    escaped << hex;

    for (string::const_iterator i = s.begin(), n = s.end(); i != n; ++i) {
        string::value_type c = (*i);

        // Keep alphanumeric and other accepted characters intact
        if (isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~') {
            escaped << c;
            continue;
        }

        // Any other characters are percent-encoded
        escaped << std::uppercase;
        escaped << '%' << std::setw(2) << int((unsigned char) c);
        escaped << std::nouppercase;
    }

    return escaped.str();
}




ostream& HttpServer::writeJQuery(ostream& html)
{
    html << "<script src='http://ajax.googleapis.com/ajax/libs/jquery/1.8.3/jquery.min.js'></script>";
    return html;

}



ostream& HttpServer::writeTableSorter(ostream& html)
{
    html << "<script src='http://tablesorter.com/__jquery.tablesorter.min.js'></script>";
    return html;
}



void HttpServer::processPost(
    const vector<string>& requestLine,
    std::iostream& s)
{
    cout << timestamp << "Received a POST." << endl;
    PostData postData(requestLine, s);
    s << "HTTP/1.1 200 OK\r\n";
    processPostRequest(postData, s);
}


PostData::PostData(const vector<string>& request, istream& s) :
    request(request)
{
    readHeaders(s);
    readContent(s);
    constructFormData();
}


void PostData::readHeaders(istream& s)
{
    // cout << "Reading POST headers." << endl;

    string line;
    while(true) {

        // Get a header line.
        getline(s, line);
        SHASTA_ASSERT(s);

        // Remove the final '\r'.
        SHASTA_ASSERT(line.size() > 0);
        SHASTA_ASSERT(line.back() == '\r');
        line.resize(line.size() - 1);

        // If empty, we have reached the end of the headers.
        if(line.size() == 0) {
            break;
        }
        // cout << "Processing POST header line:" << endl;
        // cout << line << endl;

        // Locate the ":".
        const size_t colonPosition = line.find_first_of(':');
        if(colonPosition == string::npos) {
            throw runtime_error("Missing colon in POST request.");
        }

        // Extract the keyword and the content of this header line
        // and store them.
        const string keyword = line.substr(0, colonPosition);
        const string content = line.substr(colonPosition+2);
        headers.insert(make_pair(keyword, content));
    }

}


size_t PostData::getContentLength() const
{
    // We only support POST requests that contain a content length.
    const auto it = headers.find("Content-Length");
    if(it == headers.end()) {
        throw runtime_error("POST request without content length is not supported.");
    }
    return lexical_cast<size_t>(it->second);

}

void PostData::readContent(istream& s)
{
    const size_t contentLength = getContentLength();
    // cout << "POST content length is " << contentLength << endl;

    content.resize(contentLength, '\0');
    s.read(&content.front(), content.size());

    // cout << "POST content:" << endl;
    // cout << content;
}



// Get the content boundary defined in the Content-Type header.
// The actual boundary used to separate form data in the content
// is this, prepended with two dashes.
string PostData::getBoundary() const
{
    // Locate the data in the "Content-Type" header.
    const auto it = headers.find("Content-Type");
    if(it == headers.end()) {
        throw runtime_error("POST request without content type header is not supported.");
    }
    const string& contentTypeData = it->second;

    const string target = "boundary=";
    size_t begin = contentTypeData.find(target);
    if(begin == string::npos) {
        throw runtime_error("POST request without boundary in connent type header is not supported.");
    }
    begin += target.size();
    size_t end = contentTypeData.find_first_of(' ', begin);
    if(end == string::npos) {
        end = contentTypeData.size();
    }
    return contentTypeData.substr(begin, end);

}


void PostData::constructFormData()
{
    // Get the content boundary defined in the Content-Type header.
    string boundary = getBoundary();
    // cout << "Boundary ***" << boundary << "***" << endl;

    // The actual boundary used to separate form data in the content
    // is this, prepended with two dashes.
    boundary = "--" + boundary;

    // Look for occurrences of the boundary in the content.
    vector<size_t> boundaryPositions;
    size_t startPosition = 0;
    while(startPosition < content.size()) {
        const size_t nextPosition = content.find(boundary, startPosition);
        if(nextPosition == string::npos) {
            break;
        }
        boundaryPositions.push_back(nextPosition);
        startPosition = nextPosition + boundary.size();
    }



    // Process each of the form data.
    const size_t formDataCount = boundaryPositions.size() - 1;
    // cout << "Found " << formDataCount << " form data items." << endl;
    for(size_t i=0; i<formDataCount; i++) {

        // Locate the header for this form data.
        const size_t headerBegin = boundaryPositions[i] + boundary.size();
        const size_t headerEnd = content.find("\r\n\r\n", headerBegin) + 4;

        // Locate the actual data for this form data.
        const size_t dataBegin = headerEnd;
        const size_t dataEnd = boundaryPositions[i+1] - 2; // Because of "\r\n"

        /*
        cout << "Form data item " << i << " header:" << endl;
        cout << content.substr(headerBegin, headerEnd-headerBegin);
        cout << "Form data item " << i << " data:" << endl;
        cout << content.substr(dataBegin, dataEnd-dataBegin) << endl;
        */

        // Locate the name.
        const string target = "name=\"";
        const size_t nameBegin = content.find(target, headerBegin) + target.size();
        const size_t nameEnd = content.find_first_of('"', nameBegin);
        const string name = content.substr(nameBegin, nameEnd-nameBegin);
        // cout << "Name is ***" << name << "***" << endl;

        formData.insert(make_pair(name, span<const char>(
            content.data() + dataBegin,
            content.data() + dataEnd)));

        /*
        cout << "Data length is " << formData.size() << ". Data follows:" << endl;
        for(char c: formData) {
            cout << c;
        }
        cout << endl;
        */
    }

    /*
    cout << "Received the following POST data:" << endl;
    for(const auto& p: formData) {
        cout << p.first << " of length " << p.second.size() << endl;
        for(const char c: p.second) {
            cout << c;
        }
        cout << endl;
    }
    */
}


// If the derived class wants to process POST requests, it should override this.
void HttpServer::processPostRequest(
    const PostData&,
    ostream& html)
{
    html << "POST request ignored.";
    cout << "\nPOST request ignored." << endl;
}



// Return true if the connection is a local connection
// originating from a process owned by the same
// user running the server.
// We use the lsof command to find open sockets on 127.0.0.1
// and the current port as source or destination.
bool HttpServer::isLocalConnectionSameUser(
    boost::asio::ip::tcp::iostream& s,
    uint16_t port) const
{
    // Get out process id.
    const string serverProcessId = to_string(::getpid());

    // Use lsof command to get a list of open sockets.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string fileName = tmpDirectory() + uuid + ".fa";
    if(::system(("lsof -i -n > " + fileName).c_str())) {
        filesystem::remove(fileName);
        return false;
    }



    // Read the output of the lsof command.
    ifstream file(fileName);
    string line;
    vector<string> tokens;
    const std::regex sourceAndDestinationRegex("(.+)\\-\\>(.+)");
    const std::regex colonRegex("(.+)\\:(.+)");
    string serverProcessUserName;
    string clientProcessUserName;
    const string portString = to_string(port);
    while(true) {

        // Get a line.
        getline(file, line);
        if(!file) {
            break;
        }

        // Parse it.
        boost::algorithm::split(tokens, line,
            boost::algorithm::is_any_of(" "),
            boost::token_compress_on);

        // Filter out the lines we are not interested in.
        if(tokens.size() != 10) {
            continue;
        }
        if(tokens[7] != "TCP") {
            continue;
        }
        if(tokens[9] != "(ESTABLISHED)") {
            continue;
        }

        // Extract the source and destination.
        const string& sourceAndDestination = tokens[8];
        std::smatch matches;
        const bool isMatch = std::regex_match(
            sourceAndDestination,
            matches,
            sourceAndDestinationRegex);
        if(!isMatch) {
            continue;
        }
        if(matches.size() != 3) {
            continue;
        }
        const string& source = matches[1];
        const string& destination = matches[2];

        // If the source or destination begin with "[",
        // they are IPv6 addresses and we skip them.
        // We are only interested in 127.0.0.1.
        if(source.size()==0 || source[0]=='[') {
            continue;
        }
        if(destination.size()==0 || destination[0]=='[') {
            continue;
        }

        // Parse the source and destination into IP address and port.
        if(!std::regex_match(source, matches, colonRegex)) {
            continue;
        }
        if(matches.size() != 3) {
            continue;
        }
        const string sourceIpAddress = matches[1];
        const string sourcePort = matches[2];
        if(!std::regex_match(destination, matches, colonRegex)) {
            continue;
        }
        if(matches.size() != 3) {
            continue;
        }
        const string destinationIpAddress = matches[1];
        const string destinationPort = matches[2];

        const string& processId = tokens[1];
        const string& userName = tokens[2];

        /*
        cout << line << endl;
        cout <<
            serverProcessId << " " <<
            processId << " " <<
            userName << " " <<
            sourceIpAddress << " " <<
            sourcePort << " " <<
            destinationIpAddress << " " <<
            destinationPort << " " <<
            endl;
        */

        if( sourceIpAddress == "127.0.0.1" &&
            destinationIpAddress == "127.0.0.1" &&
            sourcePort==portString &&
            processId==serverProcessId) {
            serverProcessUserName = userName;
        }
        if( sourceIpAddress == "127.0.0.1" &&
            destinationIpAddress == "127.0.0.1" &&
            destinationPort==portString) {
            clientProcessUserName = userName;
        }
    }
    filesystem::remove(fileName);

    return
        serverProcessUserName.size()>0 &&
        clientProcessUserName.size()>0 &&
        serverProcessUserName==clientProcessUserName;

    return false;    // For now.
}


#endif
