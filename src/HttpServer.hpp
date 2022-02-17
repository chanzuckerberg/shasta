#ifdef SHASTA_HTTP_SERVER

// Class http server can be used as a base class to provide simple
// http server functionality to facilitate data exploration and debugging.
// The derived class only has to override
// function processRequest.

#ifndef SHASTA_HTTP_SERVER_HPP
#define SHASTA_HTTP_SERVER_HPP

#include "span.hpp"

#include <boost/asio/ip/tcp.hpp>
#include <boost/lexical_cast.hpp>

#include "iosfwd.hpp"
#include <map>
#include <set>
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    class HttpServer;
    class PostData;
}



/*******************************************************************************

The second and third argument to HttpServer::explore control
how strict the server will be in accepting and processing connections.

If localOnly is true, only local connections are accepted
(connections from localhost, 127.0.0.1).
Remote connections (connections originating from other hosts) are not accepted.

If in addition sameUserOnly is also true, only local connections
from processes belonging to the same user are processed.
Local connections from processes belonging to different users
are initially accepted but then closed before returning any data.

If localOnly is false, sameUserOnly must also be false.

In summary, here is what the 4 combinations of localOnly and sameUserOnly do:

localOnly=false, sameUserOnly=false
All connections are accepted: local and remote, from any user.
This means that all users, not only on the computer running the server,
but also on all computers on the same local area network, can use the server.
In addition, if the computer running the server is not protected by a firewall,
everybody on the Internet can also access the server.
THIS CHOICE ALLOWS OTHER USERS TO LOOK AT YOUR DATA AND YOU SHOULD
NOT BE USED IF ACCESS TO THE DATA SHOULD BE RESTRICTED.

localOnly=false, sameUserOnly=true
This combination is invalid and results in an error.

localOnly=true, sameUserOnly=false
Only connections from the local computer are accepted, from any user.
This means that other users with access to the local computer can
access the server.
THIS CHOICE ALLOWS OTHER USERS TO LOOK AT YOUR DATA AND YOU SHOULD
NOT BE USED IF ACCESS TO THE DATA SHOULD BE RESTRICTED.

localOnly=true, sameUserOnly=true
Only connections from the local computer and originating
on a process owned by the same user running the server
are accepted. This is the only choice that limits
access to the data to the same user running the server.

*******************************************************************************/



class shasta::HttpServer {
public:

    // This function puts the server into an endless loop of processing requests.
    // See comments above for the meaning of localOnly and sameUserOnly
    void explore(uint16_t port, bool localOnly, bool sameUserOnly);

    // The destructor needs to be virtual for clean destruction of
    // the derived class.
    virtual ~HttpServer()
    {
    }

protected:

    // The derived class should override this.
    // It is passed the string of the GET request,
    // already parsed using "?=&" as separators.
    // It should write the response to the given request on the stream passed as a second argument.
    // The request is guaranteed not to be empty.
    class BrowserInformation;
    virtual void processRequest(
        const vector<string>& request,
        ostream& html,
        const BrowserInformation&) = 0;

    // If the derived class wants to process POST requests, it should override this.
    virtual void processPostRequest(
        const PostData&,
        ostream& html);


public:
    // This function can be used to get the value of a parameter.
    // If the parameter is missing, returns false and the value is not touched.
    template<class T> static bool getParameterValue(const vector<string>& request, const string& name, T& value)
    {
        for(size_t i = 0; i < request.size() - 1; i++) {
            if(request[i] == name) {
                try {
                    value = boost::lexical_cast<T>(request[i + 1]);
                } catch (...) {
                    return false;
                }
                return true;
            }
        }
        return false;
    }

    // Return all values assigned to a parameter.
    // For example, if the request has ...&a=xyz&a=uv,
    // when called with argument "a" returns a set containing "xyz" and "uv".
    static void getParameterValues(const vector<string>& request, const string& name, vector<string>& values);
    static void getParameterValues(const vector<string>& request, const string& name, std::set<string>& values);
protected:

    static ostream& writeJQuery(ostream& html);
    static ostream& writeTableSorter(ostream& html);

    // This takes care of percent encoding in the url.
    // I adapted it from boost asio example http/server3.
    // This is necessary if we want to be able to accept form data that
    // contain characters that are forbidden in an URL.
    bool urlDecode(const string& in, string& out);


    // Function to do URL encoding.
    // This is necessary when we create a hyperlink that contains special characters.
    // Adapted from here:
    // https://stackoverflow.com/questions/154536/encode-decode-urls-in-c
    static string urlEncode(const string&);


protected:

    // Class used to identify the browser that issued a request.
    class BrowserInformation {
    public:
        bool isChrome = false;
        bool isFirefox = false;
        bool isEdge = false;
        void set(const string& userAgentHeader);
    };



private:
    void processRequest(boost::asio::ip::tcp::iostream&);

    void processPost(
        const vector<string>& request,
        std::iostream&);

    void setRequestTimeout(int, boost::asio::ip::tcp::iostream&);

    // Return true if the connection is a local connection
    // originating from a process owned by the same
    // user running the server.
    bool isLocalConnectionSameUser(boost::asio::ip::tcp::iostream&, uint16_t port) const;
};



// Class describing a POST request.
class shasta::PostData {
public:

    // The request already parsed in tokens.
    vector<string> request;

    // The headers.
    std::map<string, string> headers;

    // The content data.
    string content;

    // The form data.
    // Keyed by the name.
    // The values point into the content data string above.
    std::map<string, span<const char> > formData;

    PostData(const vector<string>& request, istream&);

private:
    void readHeaders(istream&);
    size_t getContentLength() const;
    void readContent(istream&);
    void constructFormData();

    // Get the content boundary defined in the Content-Type header.
    // The actual boundary used to separate form data in the content
    // is this, prepended with two dashes.
    string getBoundary() const;


};
#endif

#endif
