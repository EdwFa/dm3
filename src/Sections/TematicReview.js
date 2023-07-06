import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'
//import "./network.css";

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';

import { variables } from './Variables.js';
import { Search } from './TematicReview/Search.js'
import { Review } from './TematicReview/Review.js'


export class TematicReview extends Component {

    constructor(props) {
        super(props);

        this.gridRef = createRef();
        this.gridAnaliseRef = createRef();
        this.state = {
            updateOr: false,
            loading: false,
            token: variables.token,
            current_page: 'search',
        }
    }
    getPage(name) {
      this.setState({current_page: name})
    }

    componentDidMount() {
        console.log('start');
    }

    render() {
        const {
            token,
            current_page,
        } = this.state;

        if (!token){
            return <Navigate push to="/login" />
        } else {
            return (
            <>
                <div className="container-fluid">
                    <div className="row">
                      <header className="navbar navbar-dark sticky-top bg-primary flex-md-nowrap p-0 shadow-sm">
                      <div className="col-md-2">
                        <div className="row g-0">
                          <div className="col-md-2">
                            <button className="mt-1 navbar-toggler collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-controls="sidebar" aria-expanded="false" aria-label="Toggle navigation" >
                              <span className="icon-bar top-bar"></span>
                              <span className="icon-bar middle-bar"></span>
                              <span className="icon-bar bottom-bar"></span>
                            </button>
                          </div>
                          <div className="col-md-10">
                            <a className="navbar-brand" href="#">SECHENOV AI-DATAMED</a>
                          </div>
                        </div>
                      </div>
                      <div className="col-md-8">
                      </div>
                      <div className="col-md-2">
                        <div className="row g-0">
                          <div className="col-md-2">
                            <button className="mt-2 navbar-toggler collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#sidebar2" aria-controls="sidebar2" aria-expanded="false" aria-label="Toggle navigation">
                              <span className="icon-bar top-bar"></span>
                              <span className="icon-bar middle-bar"></span>
                              <span className="icon-bar bottom-bar"></span>
                            </button>
                          </div>
                        </div>
                      </div>
                      </header>
                    </div>
                </div>
                <main>
                    <div>
                        <div className="container-fluid">
                            <div className="row align-items-stretch b-height">

                            <aside id="sidebar" className="col-md-2 bg-light collapse show width mb-5 shadow-sm g-0">
                                <div class="bd-example">
                                    <ul class="nav nav-pills mb-3" id="myTab" role="tablist">
                                      <li class="nav-item" role="presentation">
                                        <button class="nav-link active" id="home-tab" data-bs-toggle="tab" data-bs-target="#home" type="button" role="tab" aria-controls="home" aria-selected="false" onClick={() => this.getPage('search')}>Результаты поиска</button>
                                      </li>
                                      <li class="nav-item" role="presentation">
                                        <button class="nav-link" id="profile-tab" data-bs-toggle="tab" data-bs-target="#profile" type="button" role="tab" aria-controls="profile" aria-selected="true" onClick={() => this.getPage('review')}>Тематическое описание коллекции</button>
                                      </li>
                                      <li class="nav-item" role="presentation">
                                        <button class="nav-link" id="contact-tab" data-bs-toggle="tab" data-bs-target="#contact" type="button" role="tab" aria-controls="contact" aria-selected="false" onClick={() => this.getPage('Some')}>Схема</button>
                                      </li>
                                    </ul>
                                </div>
                            </aside>
                            <>
                                {current_page === 'search'?
                                    <Search />
                                :
                                current_page === 'review'?
                                    <Review />
                                :null
                                }
                            </>
                            </div>
                        </div>
                    </div>
                </main>
            </>
            )
        }
    }
}