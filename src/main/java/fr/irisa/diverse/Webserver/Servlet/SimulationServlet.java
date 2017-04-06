package fr.irisa.diverse.Webserver.Servlet;

import fr.irisa.diverse.Utils.Path;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.*;

/**
 * Servlet used to launch the simulation of the boussinesq model. It will run the docker image and the simulation
 * into it.
 *
 * Created by antoine on 20/03/17.
 */
public class SimulationServlet extends HttpServlet {

    public SimulationServlet() {
        // Don't do anything for now
    }

    @Override
    public void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        // Describe the script used to launch the simulation
        File script = createSimulationScript();

        // Run the docker command from a script
        Process process;
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", script.toString());

            // Redirect the input, output and error to terminal
            pb.redirectErrorStream(true);
            pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);

            process = pb.start();
            System.out.println("Simulation started");
            process.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
            response.setStatus(HttpServletResponse.SC_INTERNAL_SERVER_ERROR);
        } finally {
            // Everything went well
            script.delete();
            response.setStatus(HttpServletResponse.SC_OK);
            System.out.println("Simulation finished");
        }

    }

    // Custom method
    private File createSimulationScript() throws IOException {
        String absoluteURIToHs1d = Path.getHs1dPath();
        System.out.println(absoluteURIToHs1d);

        File script = File.createTempFile("script", null);

        Writer streamWriter = new OutputStreamWriter(new FileOutputStream(script));
        PrintWriter printWriter = new PrintWriter(streamWriter);

        printWriter.println("#!/bin/bash");
        printWriter.println("docker run --rm -v " + absoluteURIToHs1d + ":/app antoinecheronirisa/hs1d:latest /bin/bash -c ./entrypoint.sh");

        printWriter.close();

        return script;
    }

    private void displayOutputStream(OutputStream os) {
        BufferedReader reader = null;
        System.out.printf("Output Stream : " + os.toString());
    }

    private void displayErrorStream(InputStream es) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(es));
        System.out.printf("Error Stream : ");
        String line;
        while((line = reader.readLine()) != null) {
            System.out.println(line);
        }
    }
}
