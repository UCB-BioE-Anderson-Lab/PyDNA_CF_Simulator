from flask import Flask, request, jsonify
import yaml
from flask_cors import CORS
from pydna_cf_simulator.parse_CF_shorthand import parse_CF_shorthand
from pydna_cf_simulator.simulate_CF import simulate_CF

app = Flask(__name__)
CORS(app)

@app.route('/simulate', methods=['POST'])
def simulate():
    cf_shorthand = request.json['cf']
    cf = parse_CF_shorthand(cf_shorthand)
    result = simulate_CF(cf)
    return jsonify({k: vars(v) for k, v in result.items()})

@app.route('/.well-known/ai-plugin.json')
def serve_plugin_json():
    with open('ai-plugin.json', 'r') as f:
        plugin_json = yaml.safe_load(f)
    with open('docs/cf_shorthand_specification.md', 'r') as f:
        plugin_json['description_for_model'] += '\n\n' + f.read()
    return jsonify(plugin_json)

@app.route('/openapi.yaml')
def serve_openapi_yaml():
    with open('openapi.yaml', 'r') as f:
        return f.read()

if __name__ == '__main__':
    app.run(port=8234)