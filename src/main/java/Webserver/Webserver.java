package Webserver;

import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.server.handler.ContextHandler;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.handler.ResourceHandler;

/**
 *  Webserver set and configure the server
 *
 * Created by antoine on 17/03/17.
 */
public class Webserver {

    public static void main (String[] args) throws Exception {
        // Define the server
        Server server = new Server();

        // Add http connector to server
        ServerConnector http = new ServerConnector(server);
        http.setHost("localhost");
        http.setPort(8080);
        http.setIdleTimeout(30000);
        server.addConnector(http);

        // Create a new RessourceHandler that serves the static content (html, css and js files) and configures it
        ResourceHandler resourceHandler = new ResourceHandler();

        ContextHandler contextHandler1 = new ContextHandler();
        contextHandler1.setContextPath("/");
        contextHandler1.setResourceBase("src_web/");
        contextHandler1.setWelcomeFiles(new String[]{"index.html"});
        contextHandler1.setHandler(resourceHandler);

        // Add the contextHandler and the WebserverHandler to the server
        HandlerList handlers = new HandlerList();
        handlers.setHandlers(new Handler[]{ contextHandler1, new WebserverHandler()});
        server.setHandler(handlers);

        // Start the server (because it is configured :) )
        // Server.join is used to make the server join the current thread
        server.start();
        server.join();
    }
}
