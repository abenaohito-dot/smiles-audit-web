import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Base64;
import uk.ac.cam.ch.wwmm.opsin.*;

class OpsinBridge {
    static String encode(String value) {
        return Base64.getEncoder().encodeToString((value == null ? "" : value).getBytes(StandardCharsets.UTF_8));
    }
    public static void main(String[] args) throws Exception {
        BufferedReader input = new BufferedReader(new InputStreamReader(System.in, StandardCharsets.UTF_8));
        NameToStructureConfig config = new NameToStructureConfig();
        config.setWarnRatherThanFailOnUninterpretableStereochemistry(false);
        String name;
        while ((name = input.readLine()) != null) {
        OpsinResult result = NameToStructure.getInstance().parseChemicalName(name, config);
        System.out.println(result.getStatus());
        System.out.println(encode(result.getSmiles()));
        System.out.println(encode(result.getMessage()));
        }
    }
}
