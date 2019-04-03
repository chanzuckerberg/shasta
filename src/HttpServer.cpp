#ifndef SHASTA_STATIC_EXECUTABLE

// Implementation of class HttpServer - see HttpServer.hpp for more information.

#include "HttpServer.hpp"
#include "CZI_ASSERT.hpp"
#include "timestamp.hpp"
// #include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

#include <boost/algorithm/string.hpp>

#include <boost/asio/io_service.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ip/v6_only.hpp>
#include <boost/tokenizer.hpp>
using namespace boost;
using namespace asio;
using namespace ip;

#include <chrono>
#include "fstream.hpp"
#include "iostream.hpp"
#include <sstream>
#include "stdexcept.hpp"



// This function puts the server into an endless loop
// of processing requests.
// This is trhe function that the base class should call to start the server.
void HttpServer::explore(uint16_t port, bool localOnly)
{
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
        cout << "Only accepting local connections." << endl;
    } else {
        cout << "Accepting connections from all hosts." << endl;
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

          // Process the request.
          cout << timestamp << remoteEndpoint.address().to_string() << " " << flush;
          const auto t0 = std::chrono::steady_clock::now();
          processRequest(s);
          const auto t1 = std::chrono::steady_clock::now();
          const std::chrono::duration<double> t01 = t1 - t0;
          cout << timestamp << "Request satisfied in " << t01.count() << "s." << endl;
    }
}



void HttpServer::processRequest(tcp::iostream& s)
{
    // If the client is too slow sending the request, drop it.
    s.expires_from_now(boost::posix_time::seconds(1));

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
        s.expires_from_now(boost::posix_time::seconds(10000000));
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
    s.expires_from_now(boost::posix_time::seconds(86400));

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



void HttpServer::writeStyle(ostream& html)
{
    html << R"%(
<style>
    body {
        font-family: Arial;
    }
    pre {
        font-family: courier;
    }
    p, input {
        font-size: 16px;
    }
    h1, h2, h3 {
        color: DarkSlateBlue;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 2px solid MediumSlateBlue;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    td.right {
        text-align: right;
    }
    a {
        color: DarkSlateBlue;
    }
    ul.navigationMenu {
        list-style-type: none;
        margin: 0px 0px 12px 0px;
        padding: 0;
        overflow: hidden;
        background-color: #404040;
    }
    
    div.navigationButton {
        display: inline-block;
        color: white;
        text-align: center;
        padding: 14px 16px;
        text-decoration: none;
        // min-width: 120px;
    }
    
    .navigationMenuEntry:hover .navigationButton {
        background-color: black;
    }
    
    li.navigationMenuEntry {
        display: inline-block;
    }
    
    .navigationItems {
        display: none;
        position: absolute;
        background-color: DodgerBlue;
        // min-width: 120px;
        box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
        z-index: 1;
    }
    
    a.navigationItem {
        color: black;
        padding: 12px 16px;
        text-decoration: none;
        display: block;
        text-align: left;
    }
    
    .navigationItems a:hover {background-color: SteelBlue}
    
    .navigationMenuEntry:hover .navigationItems {
        display: block;
}
</style>
    )%";
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
        CZI_ASSERT(s);

        // Remove the final '\r'.
        CZI_ASSERT(line.size() > 0);
        CZI_ASSERT(line.back() == '\r');
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

        formData.insert(make_pair(name, MemoryAsContainer<const char>(
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

#endif
